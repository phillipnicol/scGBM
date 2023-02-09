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
#'
#' @return A list with components
#' \itemize{
#' \item \code{V} - A matrix containing the factor scores.
#' \item \code{U} - A matrix contianing the factor loadings
#' \item \code{D} - A vector containing the singular values (scaling factors)
#' \item \code{alpha} - A matrix containing gene-specific intercepts. The number of columns
#' is set to be equal to the number of batches (by default there are no batches so this is 1).
#' \item \code{beta} - A vector of cell-specific intercepts.
#' \item \code{I} - The number of genes
#' \item \code{J} - The number of cells.
#' \item \code{W} - The estimated mean (and by properties of the Poisson distribution, also the variance)
#' for each entry of the count matrix.
#' \item \code{obj} The value of the objective function for each iteration.
#' }
#'
#' @author Phillip B. Nicol <philnicol740@gmail.com>
gbm.sc <- function(Y,
                   M,
                   max.iter=100,
                   min.iter=0,
                   tol=10^{-4},
                   subset=NULL,
                   ncores=1,
                   infer.beta=FALSE,
                   return.W = TRUE,
                   batch=as.factor(rep(1,ncol(Y))),
                   time.by.iter = FALSE,
                   svd.free=FALSE,
                   lr=1) {
  if(!is.null(subset)) {
    out <- gbm.proj.parallel(Y,M,subsample=subset,ncores=ncores,tol=tol,
                             max.iter=max.iter)
    return(out)
  }

  I <- nrow(Y); J <- ncol(Y)
  LL <- rep(0,max.iter)
  if(time.by.iter) {
    loglik <- c()
    time <- c()
    start.time <- Sys.time()
  }

  if(!is.matrix(Y)) {
    Y <- as.matrix(Y)
  }

  if(!is.factor(batch)) {
    stop("Batch must be encoded as factor.")
  }
  batch.factor <- batch
  batch <- as.integer(batch)
  nbatch <- max(batch)

  #Precompute relevant quantities
  max.Y <- max(Y)

  #Starting estimate for alpha and W
  betas <- log(colSums(Y))
  betas <- betas - mean(betas) #Enforce the betas sum to 0
  W <- matrix(0, nrow=I, ncol=J)
  alphas <- vapply(1:nbatch, FUN.VALUE=numeric(I), function(j) {
    log(rowSums(Y[,batch==j]))-log(sum(exp(betas[batch==j])))
  })
  W <- t(t(exp(alphas[,batch]))*exp(betas))

  #Starting estimate of X
  Z <- (Y-W)/sqrt(W)
  Z[Z > sqrt(2*log(I*J/0.025))] <- sqrt(2*log(I*J/0.025))
  Z[Z < -sqrt(2*log(I*J/0.025))] <- -sqrt(2*log(I*J/0.025))
  #Z[Z > sqrt(log(I * J))] <- sqrt(log(I * J))
  #Z[Z < -sqrt(log(I * J))] <- -sqrt(log(I * J))
  LRA <-  irlba::irlba(Z,nv=M,nu=M)
  X <- LRA$u %*%(LRA$d*t(LRA$v))
  X <- sqrt(1/W)*X
  lambdav <- var(LRA$d*t(LRA$v)[,M])
  lambdau <- var(LRA$u[,M])

  print(LRA$d[1]/LRA$d[M])
  #Bound X to avoid starting too large
  #clip <- log(Y+sqrt(J*W))
  #print(clip[1,1])
  #print(clip[2,1])
  #X[X > clip] <- clip[X>clip]
  #X[X < -clip] <- -clip[X < -clip]

  #For acceleration, save previous X
  Xt <- matrix(0,nrow=I,ncol=J)

  #Prev grad
  Gt <- matrix(0,nrow=I,ncol=J)

  for(i in 1:max.iter) {
    #Reweight
    alphas <- vapply(1:nbatch, FUN.VALUE=numeric(I), function(j) {
      log(rowSums(Y[,batch==j]))-log(rowSums(exp(t(t(X)+betas))[,batch==j]))
    })
    if(infer.beta) {
      betas <- log(colSums(Y))-log(colSums(exp(alphas[,batch]+X)))
      betas <- betas - mean(betas)
    }
    W <- t(t(exp(alphas[,batch]+X))*exp(betas))

    #Prevent W from being too large (stability)
    W[W > max.Y] <- max.Y

    #Compute working variables
    Z <- X+(Y-W)/W

    ## Compute log likelihood (no normalizing constant)
    LL[i] <- sum(Y*log(W)-W)
    if(time.by.iter) {
      time <- c(time,difftime(Sys.time(),start.time,units="sec"))
      loglik <- c(loglik,LL[i])
      start.time <- Sys.time()
    }
    cat("Iteration: ", i, ". Objective=", LL[i], "\n")
    if(i > 2) {
      tau <- abs((LL[i]-LL[i-1])/LL[i-1])
      if(tau < tol & i > min.iter) {
        break
      }
    }

    ## Gradient Step
    V <- X+((i-1)/(i+2))*(X-Xt)
    Xt <- X
    w.max <- max(W)
    W <- W/w.max
    if(svd.free) {
      Q <- W*Z+(1-W)*V
      LRA <- list(u=LRA$u, v=LRA$v)
      LRA$v <- t(Q)%*%LRA$u %*% chol2inv(chol(crossprod(LRA$u)+lambdav*diag(M)))
      LRA$u <- Q %*% LRA$v %*% chol2inv(chol(crossprod(LRA$v))+lambdau*diag(M))
      X <- LRA$u %*%t(LRA$v)
    } else {
      #Adadelta
      if(i == 1) {
        Gt <- (W*(Z-X))^2
      } else{
        Gt <- 0.1*Gt + 0.9*(W*(Z-X))^2
      }

      print(mean((lr/(sqrt(0.01+Gt)))))
      LRA <- irlba::irlba(V+1/(sqrt(0.01+Gt))*W*(Z-X),nv=M)
      X <- LRA$u %*%(LRA$d*t(LRA$v))
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
  if(svd.free) {
    LRA <- irlba::irlba(X,nv=M)
  }
  out$V <- LRA$v; rownames(out$V) <- colnames(Y); colnames(out$V) <- 1:M
  out$D <- LRA$d; names(out$D) <- 1:M
  out$U <- LRA$u; rownames(out$U) <- rownames(Y); colnames(out$U) <- 1:M
  out$alpha <- alphas; rownames(out$alpha) <- rownames(Y)
  out$beta <- betas; names(out$beta) <- colnames(Y)
  out$M <- M
  out$I <- nrow(out$W); out$J <- ncol(out$W)
  out$obj <- LL[1:i]
  if(time.by.iter) {
    out$loglik <- loglik
    out$time <- cumsum(time)
  }

  out <- process.results(out)
  return(out)
}

gbm.proj.parallel <- function(Y,M,subsample=2000,min.counts=5,
                              ncores,tol=10^{-4},max.iter=max.iter) {

  J <- ncol(Y); I <- nrow(Y)
  alphas.full <- log(rowSums(Y))-log(sum(colSums(Y)))
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
  colnames(U) <- paste0("U",1:M)
  alpha <- alphas.full[ixs]

  Y <- Y[ixs,]
  print(Sys.time())
  V <- matrix(0,nrow=J,ncol=2*M)
  split.len <- 200000
  max.iter <- ceiling(J/split.len)
  for(i in 1:max.iter) {
    start <- (i-1)*split.len + 1
    if(i == max.iter) {
      stop <- J
    } else {
      stop <- i*split.len
    }
    V.sub <- mclapply(start:stop, function(j) {
      cell <- as.vector(Y[,j])
      o <- log(sum(cell))+alpha
      val <- matrix(rep(0,2*M),nrow=1)
      try({
        fit <- fastglm::fastglm(x=U,y=cell,offset=o,
                       family=poisson(),
                       method=3)
        #fit <- glm(cell~0+offset(o)+.,
        #data=U,
        #family=poisson(link="log"))
        val <- matrix(c(fit$coefficients,fit$se),
                      nrow=1)
      })
      val
    },
    mc.cores = ncores)
    V.sub <- matrix(unlist(V.sub),nrow=stop-start+1,ncol=2*M,byrow=TRUE)
    V[start:stop,] <- V.sub
  }
  out$se_V <- V[,(M+1):(2*M)]
  out$V <- V[,1:M]

  alpha <- rep(0, I)
  alpha[ixs] <- out$alpha[,1]
  alpha[-ixs] <- min(out$alpha[,1])
  out$alpha <- alphas.full
  U <- matrix(0, nrow=I,ncol=M)
  U[ixs,] <- out$U
  out$U <- U
  return(out)
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

  gbm$V <- t(diag(gbm$D) %*% t(gbm$V))

  return(gbm)
}



infer.gbm <- function(out) {
  M <- length(out$D)
  U <- out$U
  V <- out$V
  W <- out$W
  J <- nrow(V)
  tU <- t(out$U)
  se_V <- apply(W,2,function(w) {
    Mat <- tU %*% (w*out$U)
    sqrt(diag(solve(Mat)))
  })
  se_V <- t(se_V)
  return(se_V)
}


get.seu <- function(out) {
  M <- length(out$D)
  U <- out$U
  V <- out$V
  W <- out$W
  J <- nrow(V)
  tV <- t(V)
  se_U <- apply(W,1,function(w) {
    Mat <- (tV %*% (w*V))
    sqrt(diag(solve(Mat)))
  })
  se_U <- t(se_U)
  return(se_U)
}


gbm.sc2 <- function(Y,
                   M,
                   max.iter=100,
                   tol=10^{-4},
                   subset=NULL,
                   ncores=1,
                   infer.beta=FALSE,
                   return.W = TRUE,
                   batch=as.factor(rep(1,ncol(Y))),
                   time.by.iter = FALSE) {
  if(!is.null(subset)) {
    out <- gbm.proj.parallel(Y,M,subsample=subset,ncores=ncores,tol=tol,
                             max.iter=max.iter)
    return(out)
  }

  I <- nrow(Y); J <- ncol(Y)
  LL <- rep(0,max.iter)
  if(time.by.iter) {
    print("H")
    loglik <- rep(0,max.iter)
    time <- rep(0, max.iter)
    start.time <- Sys.time()
  }

  if(!is.matrix(Y)) {
    Y <- as.matrix(Y)
  }

  if(!is.factor(batch)) {
    stop("Batch must be encoded as factor.")
  }
  batch.factor <- batch
  batch <- as.integer(batch)
  nbatch <- max(batch)

  #Precompute relevant quantities
  max.Y <- max(Y)

  #Starting estimate for alpha and W
  betas <- log(colSums(Y))
  betas <- betas - mean(betas) #Enforce the betas sum to 0
  W <- matrix(0, nrow=I, ncol=J)
  alphas <- vapply(1:nbatch, FUN.VALUE=numeric(I), function(j) {
    log(rowSums(Y[,batch==j]))-log(sum(exp(betas[batch==j])))
  })
  W <- t(t(exp(alphas[,batch]))*exp(betas))

  #Starting estimate of X
  Z <- (Y-W)/sqrt(W)
  LRA <-  irlba::irlba(Z,nv=M,nu=M)
  X <- LRA$u %*%(LRA$d*t(LRA$v))
  X <- sqrt(1/W)*X

  #Bound X to avoid starting too large
  X[X > 8] <- 8
  X[X < -8] <- -8

  #For acceleration, save previous X
  Xt <- matrix(0,nrow=I,ncol=J)

  LRA$v <- t(LRA$d*t(LRA$v))
  for(i in 1:max.iter) {
    #Reweight
    alphas <- vapply(1:nbatch, FUN.VALUE=numeric(I), function(j) {
      log(rowSums(Y[,batch==j]))-log(rowSums(exp(t(t(X)+betas))[,batch==j]))
    })
    if(infer.beta) {
      betas <- log(colSums(Y))-log(colSums(exp(alphas[,batch]+X)))
      betas <- betas - mean(betas)
    }
    W <- t(t(exp(alphas[,batch]+X))*exp(betas))

    #Prevent W from being too large (stability)
    W[W > max.Y] <- max.Y

    #Compute working variables
    Z <- X+(Y-W)/W

    ## Compute log likelihood (no normalizing constant)
    LL[i] <- sum(Y*log(W)-W)
    cat("Iteration: ", i, ". Objective=", LL[i], "\n")
    if(i > 2) {
      tau <- abs((LL[i]-LL[i-1])/LL[i-1])
      if(tau < tol) {
        break
      }
    }
    if(time.by.iter) {
      time[i] <- difftime(Sys.time(),start.time,units="sec")
      loglik[i] <- LL[i]
      start.time <- Sys.time()
    }

    ## Gradient Step
    V <- X+((i-1)/(i+2))*(X-Xt)
    Xt <- X
    w.max <- max(W)
    W <- W/w.max
    Q <- W*Z+(1-W)*V
    LRA <- list(u=LRA$u)
    LRA$v <- t(Q)%*%LRA$u %*% chol2inv(chol(crossprod(LRA$u)))
    LRA$u <- Q %*% LRA$v %*% chol2inv(chol(crossprod(LRA$v)))
    X <- LRA$u %*%t(LRA$v)

    if(i == max.iter) {
      warning("Maximum number of iterations reached (increase max.iter).
              Possible non-convergence.")
    }
  }

  LRA <- irlba::irlba(X,nv=M)
  out <- list()
  if(return.W) {
    out$W <- t(t(exp(alphas[,batch]+X))*exp(betas))
  }
  out$V <- LRA$v; rownames(out$V) <- colnames(Y); colnames(out$V) <- 1:M
  out$D <- LRA$d; names(out$D) <- 1:M
  out$U <- LRA$u; rownames(out$U) <- rownames(Y); colnames(out$U) <- 1:M
  out$alpha <- alphas; rownames(out$alpha) <- rownames(Y)
  out$beta <- betas; names(out$beta) <- colnames(Y)
  out$M <- M
  out$I <- nrow(out$W); out$J <- ncol(out$W)
  out$obj <- LL[1:i]
  if(time.by.iter) {
    out$loglik <- loglik
    out$time <- cumsum(time)
  }

  out <- process.results(out)
  return(out)
}
