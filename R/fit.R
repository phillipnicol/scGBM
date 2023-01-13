

gbm.sc <- function(Y,
                   M,
                   max.iter=100,
                   tol=10^{-4},
                   subset=NULL,
                   ncores=1,
                   infer.beta=FALSE,
                   return.W = TRUE,
                   batch=as.factor(rep(1,ncol(Y)))) {
  if(!is.null(subset)) {
    out <- gbm.proj.parallel(Y,M,subsample=subset,ncores=ncores)
    return(out)
  }

  I <- nrow(Y); J <- ncol(Y)
  LL <- rep(0,max.iter)

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

    ## Gradient Step
    V <- X+((i-1)/(i+2))*(X-Xt)
    Xt <- X
    w.max <- max(W)
    W <- W/w.max
    LRA <- irlba::irlba(W*Z+(1-W)*V,nv=M)
    X <- LRA$u %*%(LRA$d*t(LRA$v))
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
  out$obj <- LL[1:i]

  out <- process.results(out)
  return(out)
}

gbm.proj.parallel <- function(Y,M,subsample=2000,min.counts=5,
                              ncores) {
  J <- ncol(Y)
  if(length(subsample)==1) {
    jxs <- sample(1:J,size=subsample,replace=FALSE)
    Y.sub <- Y[,jxs]
  } else {
    Y.sub <- Y[,subsample]
  }
  Y.sub <- as.matrix(Y.sub)
  ixs <- which(rowSums(Y.sub) > 5)
  Y.sub <- Y.sub[ixs,]
  out <- gbm.sc(Y.sub,M=M)

  U <- out$U
  #U <- as.data.frame(U)
  colnames(U) <- paste0("U",1:M)
  alpha <- out$alpha

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

