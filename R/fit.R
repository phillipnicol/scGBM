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
                   tol=10^{-4},
                   subset=NULL,
                   ncores=1,
                   infer.beta=FALSE,
                   return.W = TRUE,
                   batch=as.factor(rep(1,ncol(Y))),
                   time.by.iter = FALSE,
                   lr=1,
                   min.iter=30) {
  if(!is.null(subset)) {
    out <- gbm.proj.parallel(Y,M,subsample=subset,ncores=ncores,tol=tol,
                             max.iter=max.iter)
    return(out)
  }


  I <- nrow(Y); J <- ncol(Y)
  LL <- rep(0,max.iter)
  loglik <- c()
  if(time.by.iter) {
    time <- c()
    start.time <- Sys.time()
  }

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
  #betas <- betas - mean(betas) Enforce the betas sum to 0
  W <- matrix(0, nrow=I, ncol=J)
  log.rsy <- matrix(0,nrow=nbatch,ncol=I)
  for(j in 1:nbatch) {
    log.rsy[j,] <- log(rowSums(Y[,batch==j]))
  }
  alphas <- vapply(1:nbatch, FUN.VALUE=numeric(I), function(j) {
    log.rsy[j,]-log(sum(exp(betas[batch==j])))
  })
  betas <- betas + mean(alphas)
  alphas <- alphas - mean(alphas)
  W <- exp(sweep(alphas[,batch], 2, betas, "+"))

  #Starting estimate of X
  Z <- (Y-W)/sqrt(W)
  c <- sqrt(2*log(I*J/0.025))
  #Z[Z > c] <- c
  #Z[Z < -c] <- -c
  LRA <-  irlba::irlba(Z,nv=M,nu=M)
  X <- LRA$u %*%(LRA$d*t(LRA$v))
  X <- sqrt(1/W)*X

  print(sum(abs(X) > 8))
  X[X > 8] <- 8
  X[X < -8] <- -8

  #For acceleration, save previous X
  Xt <- matrix(0,nrow=I,ncol=J)

  for(i in 1:max.iter) {
    #Reweight
    alphas <- vapply(1:nbatch, FUN.VALUE=numeric(I), function(j) {
      #sweep(X[,batch==j],2,betas[batch==j],"+")
      log.rsy[j,]-log(rowSums(exp(sweep(X[,batch==j],2,betas[batch==j],"+"))))
    })
    betas <- betas + mean(alphas)
    alphas <- alphas - mean(alphas)
    if(infer.beta) {
      betas <- log(colSums(Y))-log(colSums(exp(alphas[,batch]+X)))
    }
    W <- exp(sweep(alphas[,batch]+X, 2, betas, "+"))

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

      if(LL[i] <= (LL[i-1]+0.1)) {
        lr <- max(lr/2, 1)
        X <- Xt
        next
      } else {
        lr <- lr*(1.05)
        #lr <- lr
      }
    }

    loglik <- c(loglik,LL[i])
    cat("Iteration: ", i, ". Objective=", LL[i], "lr: ", lr, "\n")
    if(time.by.iter) {
      time <- c(time,difftime(Sys.time(),start.time,units="sec"))
      start.time <- Sys.time()
    }

    ## Gradient Step
    V <- X+((i-1)/(i+2))*(X-Xt)
    Xt <- X
    w.max <- max(W)

    LRA <- irlba::irlba(V+(lr/w.max)*(Y-W),nv=M)
    X <- LRA$u %*%(LRA$d*t(LRA$v))

    if(i == max.iter) {
      warning("Maximum number of iterations reached (increase max.iter).
              Possible non-convergence.")
    }
  }

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
  out$loglik <- loglik
  if(time.by.iter) {
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
  U <- cbind(rep(1,length(ixs)),U)
  colnames(U) <- c("intercept",paste0("U",1:M))
  #alpha <- alphas.full
  alpha <- out$alpha[,1]
  #alpha <- alphas.full[ixs]

  Y <- Y[ixs,]
  print(Sys.time())
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
    V.sub <- mclapply(start:stop, function(j) {
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
  out$se_V <- V[,(M+2):(2*M+1)]
  out$V <- V[,2:(M+1)]
  out$beta <- V[,1]

  alpha <- rep(0, I)
  alpha[ixs] <- out$alpha[,1]
  alpha[-ixs] <- min(out$alpha[,1])
  out$alpha <- alpha
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
gbm.scSparse <- function(Y,
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
                   svd.free=FALSE) {
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

  if(!is.factor(batch)) {
    stop("Batch must be encoded as factor.")
  }
  batch.factor <- batch
  batch <- as.integer(batch)
  nbatch <- max(batch)

  #Precompute relevant quantities
  max.Y <- max(Y)
  nz <- which(Y != 0, arr.ind=TRUE) #Find non-zeros
  log.rsY <- log(rowSums(Y))

  #Starting estimate for alpha and W
  betas <- log(colSums(Y))
  betas <- betas - mean(betas) #Enforce the betas sum to 0

  alphas <- vapply(1:nbatch, FUN.VALUE=numeric(I), function(j) {
    log(rowSums(Y[,batch==j]))-log(sum(exp(betas[batch==j])))
  })
  alphas <- as.vector(alphas)

  #W <- Y #In sparse format
  #nz <- which(Y != 0, arr.ind=TRUE) #Find non-zeros
  #w <- alphas[nz[,1]]+betas[nz[,2]]
  #W@x <- exp(w)

  #S <- Y
  #S@x <- (Y@x - W@x)/sqrt(W@x)
  #c <- sqrt(2*log(I*J/0.025))
  #S@x[S@x > c] <- c
  #S@x[S@x < -c] <- -c

  #lra <- irlba::irlba(S,nv=M)

  #U <- as.matrix(lra$u)
  #V <- as.matrix(t(lra$d*t(lra$v)))

  W <- outer(exp(alphas),exp(betas))

  #Starting estimate of X
  Z <- (Y-W)/sqrt(W)
  Z[Z > sqrt(2*log(I*J/0.025))] <- sqrt(2*log(I*J/0.025))
  Z[Z < -sqrt(2*log(I*J/0.025))] <- -sqrt(2*log(I*J/0.025))
  #Z[Z > sqrt(log(I * J))] <- sqrt(log(I * J))
  #Z[Z < -sqrt(log(I * J))] <- -sqrt(log(I * J))
  lra <-  irlba::irlba(Z,nv=M,nu=M)
  U <- lra$u
  V <- t(lra$d*t(lra$v))

  U <- exp(alphas)^{-1/2}*U
  V <- exp(betas)^{-1/2}*V

  W <- Y
  S <- Y
  for(i in 1:max.iter) {
    #Vt1J <- colSums(V)
    #UVt1J <- as.vector(U%*%Vt1J)

    #alphas <- log.rsY - log(J+UVt1J)
    #alphas <- vapply(UVt1J, FUN.VALUE=numeric(1), function(ix) {
    #  -sum(log(exp(ix+betas)))
    #})
    #alphas <- log.rsY+alphas
    #alphas <- log.rsY - log(rowSums(exp(sweep(X,2,betas,"+"))))
    #alphas <- log(rowSums(Y)) - log(rowSums(exp(alphas)^{-1}*W))

    w <- alphas[nz[,1]]+betas[nz[,2]]+rowSums(U[nz[,1],]*V[nz[,2],])

    W@x <- exp(w)
    W@x[W@x > max.Y] <- max.Y

    #Update intercepts

    LL <- sum(Y@x*log(W@x)-W@x)
    print(LL)

    S@x <- (max(W@x))^{-1}*(Y@x - W@x)

    HU <- U%*%chol2inv(chol(crossprod(U)))
    V <- V + t(S)%*%HU

    HV <- V%*%chol2inv(chol(crossprod(V)))
    U <- U + S%*%HV
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
gbm.sc2 <- function(Y,
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
                   svd.free=FALSE) {
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

  U <- LRA$u
  V <- t(LRA$d * t(LRA$v))
  U <- exp(as.vector(alphas))^{-1/2}*U
  V <- exp(betas)^{-1/2}*V

  ixs <- which(abs(Y - W) > 0.01*max(W) ,arr.ind=TRUE)
  S <- sparseMatrix(i=ixs[,1],j=ixs[,2],x=0,dims=c(I,J))
  S@x <- (Y-W)[I*(ixs[,2]-1)+ixs[,1]]
  for(i in 1:max.iter) {
    print(i)
    #Reweight
    jxs <- sample(1:J, size=1000,replace=FALSE)
    V.sub <- V[jxs,]
    X.sub <- U %*% t(V.sub)
    alphas <- log(rowSums(Y[,jxs]))-log(rowSums(exp(sweep(X.sub,2,betas[jxs],"+"))))
    alphas[is.infinite(alphas)] <- min(alphas[!is.infinite(alphas)])
    if(infer.beta) {
      betas <- log(colSums(Y))-log(colSums(exp(alphas[,batch]+X)))
      betas <- betas - mean(betas)
    }

    w <- alphas[ixs[,1]]+betas[ixs[,2]]+rowSums(U[ixs[,1],]*V[ixs[,2],])
    W <- S
    W@x <- exp(w)
    W@x[exp(w) > max.Y] <- max.Y

    if(time.by.iter) {
      time <- c(time,difftime(Sys.time(),start.time,units="sec"))
      loglik <- c(loglik,LL[i])
      start.time <- Sys.time()
    }

    w.max <- max(W)
    #W <- W/w.max
    if(svd.free) {
      ixs2 <- which(abs(Y-W) > 0.01*w.max,arr.ind=TRUE)
      ixs <- as.data.frame(ixs); ixs2 <- as.data.frame(ixs2)
      ixs <- inner_join(ixs,ixs2)

      if(nrow(ixs) == 0) {
        break
      }
      print(nrow(ixs)/(I*J))

      ix.flat <- I*(ixs[,2]-1)+ixs[,1]
      S <- sparseMatrix(i=ixs[,1],j=ixs[,2],x=0,dims=c(I,J))
      S@x <- (w.max)^{-1}*(Y-W)[ix.flat]

      HU <- U%*%chol2inv(chol(crossprod(U)))
      V <- V + t(S)%*%HU

      HV <- V%*%chol2inv(chol(crossprod(V)))
      U <- U + S%*%HV

      X <- U %*% t(V)
    } else {
      #Adadelta
      #if(i == 1) {
      #  Gt <- (W*(Z-X))^2
      #} else{
      #  Gt <- 0.1*Gt + 0.9*(W*(Z-X))^2
      #}
      #print(mean((lr/(sqrt(1e-8+Gt)))))
      #LRA <- irlba::irlba(V+0.001/(sqrt(1e-8+Gt))*W*(Z-X),nv=M)
      print(sum(abs(Y-W) < 0.01)/(I*J))
      LRA <- irlba::irlba(V+W*(Z-X),nv=M)
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

