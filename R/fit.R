#' @export
#'
#' @title scGBM: Model-based dimensionality reduction for single-cell
#' RNA-seq with generalized bilinear models
#'
#' @description Fit a Poisson bilinear model to the single-cell
#' count matrix.
#'
#' @param Y A matrix of UMI counts with genes on rows and cells on columns.
#' @param M The number of latent factors to estimate
#' @param max.iter The maximum number of iterations
#' @param tol The tolerance for convergence (relative difference in objective function)
#' @param subset If NULL, use the entire data to fit the model. If an integer, take a
#' random subsample of cells to fit the loadings. Then project the remaining
#' cells onto those loadings. If a non-random sample is desired, subset can
#' be an integer vector corresponding to the columns that should be used to
#' fit the loadings.
#' @param ncores If subset is not NULL, then ncores specifies the number of cores to use
#' during the projection step.
#' @param infer.beta If FALSE, beta is set to the log of the number of counts per cell.
#' @param return.W If TRUE, the matrix of weights (which is equal to the estimated mean) is returned. This should
#' be FALSE for large datasets since it requires a lot of memory.
#' @param batch An optional factor containing the assignment of cells to known batches.
#' @param time.by.iter If TRUE, the elapsed time (in seconds) is given at each iteration of the algorithm.
#' @param min.iter The minimum number of iterations.
#'
#' @return A list with components
#' \itemize{
#' \item \code{scores} - A matrix containing the factor scores (weight of each cell for each factor).
#' \item \code{loadings} - A matrix (genes x factors) containing the factor loadings (weight of each cell for each factor).
#' \item \code{V} - A matrix containing the unscaled factor scores.
#' \item \code{U} - Factor loadings. This matrix is identical to \code{loadings} and is included
#' for reverse compatibility.
#' \item \code{D} - A vector containing the singular values (scaling factors)
#' \item \code{alpha} - A matrix containing gene-specific intercepts. The number of columns
#' is set to be equal to the number of batches (by default there are no batches so this is 1).
#' \item \code{beta} - A vector of cell-specific intercepts.
#' \item \code{I} - The number of genes.
#' \item \code{J} - The number of cells.
#' \item \code{W} - The estimated mean (and by properties of the Poisson distribution, also the variance)
#' for each entry of the count matrix.
#' \item \code{obj} - The value of the objective function for each iteration.
#' }
#'
#' @details scGBM fits the following model to the \eqn{I \times J} count matrix \eqn{Y}:
#' \deqn{Y \sim \text{Poisson}(\mu)}
#' \deqn{\log \mu = \alpha 1_J^T + 1_I \beta^T + U \Sigma V^T}
#' Here \eqn{\alpha} are gene-specific intercepts, \eqn{\beta_j} are cell-specific intercepts, \eqn{U} are the \eqn{I \times M} factor loadings
#' (linear combinations of genes that define factors), and \eqn{V \Sigma} are the \eqn{J \times M} scores for each cell.
#'
#' @author Phillip B. Nicol <philnicol740@gmail.com>
gbm.sc <- function(Y,
                   M,
                   max.iter=100,
                   tol=10^{-4},
                   svd_method="rsvd",
                   subset=NULL,
                   ncores=1,
                   infer.beta=FALSE,
                   return.W = TRUE,
                   batch=as.factor(rep(1,ncol(Y))),
                   time.by.iter = FALSE,
                   min.iter=30) {

  #Check validity of input
  gbm.sc.check.valid.input(as.list(environment()))

  if(!is.null(subset)) {
    out <- gbm.proj.parallel(Y,M,subsample=subset,ncores=ncores,tol=tol,
                             max.iter=max.iter)
    ##Message for users of new version about scores
    message("For users of newer versions (1.0.1+): the `scores` matrix now contains factor scores.")
    return(out)
  }

  Q.tau <- 1
  Q.f <- NULL
  Q <- NULL

  I <- nrow(Y); J <- ncol(Y)
  LL <- rep(0,max.iter)
  loglik <- c()
  if(time.by.iter) {
    time <- c()
    start.time <- Sys.time()
  }

  lr <- 1 #Default learning rate is set to 1

  if(!is.factor(batch)) {
    stop("Batch must be encoded as factor.")
  }
  batch.factor <- batch
  batch <- as.integer(batch)
  nbatch <- max(batch)

  #Precompute relevant quantities
  max.Y <- max(Y)
  nz <- which(Y != 0)

  #Starting estimate for alpha and W
  betas <- log(colSums(Y))

  W <- matrix(0, nrow=I, ncol=J)
  log.rsy <- matrix(0,nrow=nbatch,ncol=I)
  for(j in 1:nbatch) {
    log.rsy[j,] <- log(rowSums(Y[,batch==j]))
  }
  alphas <- vapply(1:nbatch, FUN.VALUE=numeric(I), function(j) {
    log.rsy[j,]-log(sum(exp(betas[batch==j])))
  })
  betas <- betas + mean(alphas)
  alphas <- alphas - mean(alphas) #Ensure alphas sum to 0
  W <- exp(sweep(alphas[,batch], 2, betas, "+"))

  #Starting estimate of X
  Z <- (Y-W)/sqrt(W)
  LRA <-  irlba::irlba(Z,nv=M,nu=M)
  X <- LRA$u %*%(LRA$d*t(LRA$v))
  X <- sqrt(1/W)*X

  X[X > 8] <- 8 #Clipping
  X[X < -8] <- -8

  #For acceleration, save previous X
  Xt <- matrix(0,nrow=I,ncol=J)

  for(i in 1:max.iter) {
    #Reweight
    Mu <- exp(Rfast::eachrow(X,betas,oper="+"))
    alphas <- vapply(1:nbatch, FUN.VALUE=numeric(I), function(j) {
      #sweep(X[,batch==j],2,betas[batch==j],"+")
      log.rsy[j,]-log(rowSums(Mu[,batch==j]))
    })
    betas <- betas + mean(alphas)
    alphas <- alphas - mean(alphas)
    if(infer.beta) {
      betas <- log(colSums(Y))-log(colSums(exp(alphas[,batch]+X)))
    }
    W <- exp(alphas[,batch])*Mu
    #W <- exp(sweep(alphas[,batch]+X, 2, betas, "+"))

    #Prevent W from being too large (stability)
    W[W > max.Y] <- max.Y

    #Compute working variables
    #Z <- X+(Y-W)/W

    ## Compute log likelihood (no normalizing constant)
    LL[i] <- sum(Y[nz]*log(W[nz]))-sum(W)
    if(is.na(LL[i]) | is.infinite(LL[i])) {
      X <- Xt
      lr <- lr/2
      i <- i - 1
      next
    }
    if(i >= 3) {
      tau <- abs((LL[i]-LL[i-2])/LL[i])
      if(tau < tol & lr <= 1.06 & i >= min.iter) {
        break
      }

      if(LL[i] < LL[i-1]) {
        lr <- max(lr/2, 1)
        X <- Xt
        next
      } else {
        lr <- lr*(1.05)
        #lr <- lr
      }
    }

    loglik <- c(loglik,LL[i])
    cat("Iteration: ", i, ". Objective=", LL[i], "\n")
    if(time.by.iter) {
      time <- c(time,difftime(Sys.time(),start.time,units="sec"))
      start.time <- Sys.time()
    }

    if(svd_method == "rsvd") {
      pgd <- pgd_rsvd(X, Xt, i, lr, W, Y, Q, Q.tau, Q.f, M)
      X <- pgd$X
      Xt <- pgd$Xt
      Q <- pgd$Q
      Q.tau <- pgd$Q.tau
      Q.f <- pgd$Q.f
    } else if(svd_method == "irlba") {
      pgd <- pgd_irlba(X, Xt, i, lr, W, Y, M)
      X <- pgd$X
      Xt <- pgd$Xt
    }

    if(i == max.iter) {
      warning("Maximum number of iterations reached (increase max.iter).
              Possible non-convergence.")
    }
  }

  out <- list()
  if(return.W) {
    out$W <- t(t(exp(alphas[,batch]+X))*exp(betas))
  }
  LRA <- pgd$LRA
  out$V <- LRA$v; rownames(out$V) <- colnames(Y); colnames(out$V) <- 1:M
  out$D <- LRA$d; names(out$D) <- 1:M
  out$U <- LRA$u; rownames(out$U) <- rownames(Y); colnames(out$U) <- 1:M
  out$loadings <- out$U
  out$alpha <- alphas; rownames(out$alpha) <- rownames(Y)
  out$beta <- betas; names(out$beta) <- colnames(Y)
  out$M <- M
  out$I <- nrow(out$W); out$J <- ncol(out$W)
  out$obj <- loglik
  if(time.by.iter) {
    out$time <- cumsum(time)
  }

  out <- process.results(out)

  ##Message for users of new version about scores
  message("For users of newer versions (1.0.1+): the `scores` matrix now contains factor scores, the `V` matrix is UNSCALED scores.")
  return(out)
}

gbm.proj.parallel <- function(Y,M,subsample=2000,min.counts=5,
                              ncores,tol=10^{-4},max.iter=max.iter) {

  J <- ncol(Y); I <- nrow(Y)
  alphas.full <- log(rowSums(Y))
  if(length(subsample)==1) {
    jxs <- sample(1:J,size=subsample,replace=FALSE)
    Y.sub <- Y[,jxs]
  } else {
    Y.sub <- Y[,subsample]
  }
  Y.sub <- as.matrix(Y.sub)
  ixs <- which(rowSums(Y.sub) > 5)
  Y.sub <- Y.sub[ixs,]
  out <- gbm.sc(Y.sub,M=M,tol=tol,max.iter=max.iter)

  U <- out$U
  #U <- as.data.frame(U)
  U <- cbind(rep(1,length(ixs)),U)
  colnames(U) <- c("intercept",paste0("U",1:M))
  #alpha <- alphas.full
  alpha <- out$alpha[,1]
  #alpha <- alphas.full[ixs]

  Y <- Y[ixs,]
  V <- matrix(0,nrow=J,ncol=2*M+1)
  split.len <- 200000
  max.iter <- ceiling(J/split.len)
  for(i in 1:max.iter) {
    start <- (i-1)*split.len + 1
    if(i == max.iter) {
      stop <- J
    } else {
      stop <- i*split.len
    }
    V.sub <- parallel::mclapply(start:stop, function(j) {
      cell <- as.vector(Y[,j])
      o <- alpha#+log(sum(cell))
      val <- matrix(rep(0,2*M+1),nrow=1)
      try({
        fit <- fastglm::fastglm(x=U,y=cell,offset=o,
                       family=poisson(),
                       method=3)
        #fit <- glm(cell~0+offset(o)+.,
        #data=U,
        #family=poisson(link="log"))
        val <- matrix(c(fit$coefficients,fit$se[-1]),
                      nrow=1)
      })
      val
    },
    mc.cores = ncores)
    V.sub <- matrix(unlist(V.sub),nrow=stop-start+1,ncol=2*M+1,byrow=TRUE)
    V[start:stop,] <- V.sub
  }
  out$se_scores <- V[,(M+2):(2*M+1)]
  out$scores <- V[,2:(M+1)]
  out$beta <- V[,1]

  alpha <- rep(0, I)
  alpha[ixs] <- out$alpha[,1]
  alpha[-ixs] <- alphas.full[-ixs] - log(sum(exp(out$beta)))
  out$beta <- out$beta + mean(alpha)
  alpha <- alpha - mean(alpha)
  out$alpha <- alpha
  U <- matrix(0, nrow=I,ncol=M)
  U[ixs,] <- out$U
  out$U <- U
  out$V <- NULL #In projection method, don't return V.
  out$loadings <- loadings

  return(out)
}

gbm.sc.check.valid.input <- function(my.args) {
  ### Check count matrix for NA or negative values
  if(sum(is.na(my.args$Y))) {
    stop("Count matrix Y cannot contain NA values.")
  }
  else if(min(my.args$Y) < 0) {
    stop("Count matrix Y cannot contain negative values.")
  }
  ### Check M
  if(my.args$M < 1) {
    stop("M must be >= 1")
  }
}

process.results <- function(gbm) {
  #Enforce identifiability in U
  M <- gbm$M
  for(m in 1:M) {
    if(gbm$U[1,m] < 0) {
      gbm$U[,m] <- -1*gbm$U[,m]
      gbm$V[,m] <- -1*gbm$V[,m]
    }
  }
  gbm$scores <- t(gbm$D*t(gbm$V))

  return(gbm)
}

pgd_irlba <- function(X,Xt,i,lr,W,Y,M) {
  out <- list()

  ## Gradient Step
  V <- X+((i-1)/(i+2))*(X-Xt)
  out$Xt <- X
  w.max <- max(W)

  LRA <- irlba::irlba(V+(lr/w.max)*(Y-W),nv=M)
  out$X <- LRA$u %*% (LRA$d*t(LRA$v))
  out$LRA <- LRA

  return(out)
}

pgd_rsvd <- function(X, Xt, i, lr, W, Y, Q, Q.tau, Q.f, M) {
  out <- list()

  ## Gradient Step
  V <- X+((i-1)/(i+2))*(X-Xt)
  out$Xt <- X
  w.max <- max(W)

  if(i >= 2) {
    if(Q.tau < 0.001) {
      B <- crossprod(Q.f, V+(lr/w.max)*(Y-W))
      LRA <- svd(B)
      LRA$u <- Q.f %*% LRA$u
      LRA$u <- LRA$u[,1:M]
      LRA$v <- LRA$v[,1:M]
      LRA$d <- LRA$d[1:M]
    } else {
      LRA <- irlba::svdr(V+(lr/w.max)*(Y-W), Q=Q, k=M,return.Q=TRUE, it=2L,
                         extra=30)
      Q.new <- LRA$Q
      Q.tau <- mean(abs(Q-Q.new))
      print(Q.tau)
      Q <- Q.new
      Q.f <- qr.Q(qr((V+(lr/w.max)*(Y-W)) %*% Q.new))
    }
  } else {
    LRA <- irlba::svdr(V+(lr/w.max)*(Y-W),k=M, return.Q=TRUE, it=3L,
                       extra=30)
    Q <- LRA$Q
  }

  out$X <- LRA$u %*% (LRA$d*t(LRA$v))
  out$LRA <- LRA
  out$Q <- Q
  out$Q.f <- Q.f
  out$Q.tau <- Q.tau
  return(out)
}
