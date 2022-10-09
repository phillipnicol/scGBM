

gbm.sc <- function(Y,
                   M,
                   max.iter=100,
                   tol=10^{-4},
                   subset=NULL,
                   ncores=1,
                   infer.beta=FALSE,
                   big.matrix=FALSE) {
  if(!is.null(subset)) {
    if(ncores==1) {
      out <- gbm.projection(Y,M,subsample=subset)
    } else{
      out <- gbm.proj.parallel(Y,M,subsample=subset,ncores=ncores)
    }
    return(out)
  }

  I <- nrow(Y); J <- ncol(Y)
  LL <- rep(0,max.iter)

  if(!is.matrix(Y)) {
    Y <- as.matrix(Y)
  }

  #Precompute relevant quantities
  max.Y <- max(Y)


  #Starting estimate for alpha and W
  betas <- log(colSums(Y))
  alphas <- log(rowSums(Y))-log(sum(exp(betas)))
  W <- outer(exp(alphas),exp(betas))

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
    alphas <- log(rowSums(Y))-log(rowSums(exp(t(t(X)+betas))))
    if(infer.beta) {
      betas <- log(colSums(Y))-log(colSums(exp(alphas+X)))
    }
    W <- outer(exp(alphas),exp(betas))*exp(X)

    #Prevent W from being too large (stability)
    W[W > max.Y] <- max(Y)

    #Compute working variables
    Z <- X+(Y-W)/W

    ## Compute log likelihood
    LL[i] <- sum(dpois(Y,lambda=W,log=TRUE))
    cat("Iteration: ", i, ". LogLik=", LL[i], "\n")
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
  out$W <- outer(exp(alphas),exp(betas))*exp(X)
  out$V <- LRA$v
  out$D <- LRA$d
  out$U <- LRA$u
  out$alpha <- alphas
  out$beta <- betas
  out$M <- M
  out$I <- nrow(out$W); out$J <- ncol(out$W)
  out$LL <- LL

  out <- process.results(out)
  return(out)
}


gbm.projection <- function(Y,M,subsample=2000,min.counts=5,ncores) {
  J <- ncol(Y)
  jxs <- sample(1:J,size=subsample,replace=FALSE)
  Y.sub <- Y[,jxs]
  Y.sub <- as.matrix(Y.sub)
  ixs <- which(rowSums(Y.sub) > 5)
  Y.sub <- Y.sub[ixs,]
  out <- gbm.sc(Y.sub,M=M)

  V.subsamp <- t(diag(out$D) %*% t(out$V))

  U <- out$U
  U <- as.data.frame(U)
  colnames(U) <- paste0("U",1:M)
  alpha <- out$alpha

  Y <- Y[ixs,]
  #Now get V
  V <- apply(Y,2,function(cell) {
    o <- log(sum(cell))+alpha
    fit <- glm(cell~0+offset(o)+.,
               data=U,
               family=poisson(link="log"))
    sumfit <- summary(fit)
    c(coefficients(fit),sumfit$coefficients[,2])
  })
  out$se_V <- V[,(M+1):(2*M)]
  out$V <- V[,1:M]

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
  V <- mclapply(1:J,function(j) {
    cell <- as.vector(Y[,j])
    o <- log(sum(cell))+alpha
    val <- matrix(rep(0,2*M),nrow=1)
    try({
      fit <- fastglm(x=U,y=cell,offset=o,
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
  V <- matrix(unlist(V),nrow=J,ncol=2*M,byrow=TRUE)

  out$se_V <- V[,(M+1):(2*M)]
  out$V <- V[,1:M]

  return(out)
}


gbm.projection <- function(Y,M,subsample=2000,min.counts=5,ncores) {
  J <- ncol(Y)
  jxs <- sample(1:J,size=subsample,replace=FALSE)
  Y.sub <- Y[,jxs]
  Y.sub <- as.matrix(Y.sub)
  ixs <- which(rowSums(Y.sub) > 5)
  Y.sub <- Y.sub[ixs,]
  out <- gbm.sc(Y.sub,M=M)

  V.subsamp <- t(diag(out$D) %*% t(out$V))

  U <- out$U
  U <- as.data.frame(U)
  colnames(U) <- paste0("U",1:M)
  alpha <- out$alpha

  Y <- Y[ixs,]
  #Now get V
  V <- apply(Y,2,function(cell) {
    o <- log(sum(cell))+alpha
    fit <- glm(cell~0+offset(o)+.,
               data=U,
               family=poisson(link="log"))
    sumfit <- summary(fit)
    c(coefficients(fit),sumfit$coefficients[,2])
  })
  out$se_V <- V[,(M+1):(2*M)]
  out$V <- V[,1:M]

  return(out)
}





gbm.sc.batch <- function(Y,M,max.iter=100,tol=10^{-4},
                   offset=NULL,subset=NULL,ncores=1,
                   batch=rep(1,ncol(Y))) {
  if(!is.null(subset)) {
    if(ncores==1) {
      out <- gbm.projection(Y,M,subsample=subset)
    } else{
      out <- gbm.proj.parallel(Y,M,subsample=subset,ncores=ncores)
    }
    return(out)
  }


  if(!is.matrix(Y)) {
    Y <- as.matrix(Y)
  }

  nbatch <- max(batch)
  if(is.null(offset)) {
    offset <- colSums(Y)
  }

  I <- nrow(Y); J <- ncol(Y)
  LL <- rep(0,max.iter)

  #Starting estimate for mu
  mu <- vapply(1:nbatch,FUN.VALUE=numeric(I),function(i) {
    rowSums(Y[,batch==i])/sum(offset[batch==i])
  })
  cut <- apply(mu,1,function(x) {
    if(min(x) == 0) {
      1
    } else{
      0
    }
  })
  cut <- which(cut==1)
  Y <- Y[-cut,]
  mu <- mu[-cut,]
  I <- nrow(Y); J <- ncol(Y)
  alphas <- log(mu)
  Offset <- matrix(offset,nrow=I,ncol=J,byrow=TRUE)
  W <- matrix(0,nrow=I,ncol=J)
  for(i in 1:nbatch) {
    W[,batch==i] <- matrix(mu[,i], I, sum(batch==i))
  }
  W <- W*Offset

  Z <- (Y-W)/sqrt(W)
  LRA <-  irlba::irlba(Z,nv=M,nu=M)
  X <- LRA$u %*% diag(LRA$d) %*% t(LRA$v)
  X <- sqrt(1/W)*X
  LRA <- irlba::irlba(X,nv=M,nu=M)

  #Bound X to avoid starting too large
  X[X > 8] <- 8
  X[X < -8] <- -8

  Xt <- matrix(0,nrow=I,ncol=J)

  last.change <- 1

  for(i in 1:max.iter) {
    #Reweight
    Offset.new <- Offset*exp(X)
    mu <- rowSums(Y)/rowSums(Offset.new)
    mu <- vapply(1:nbatch,FUN.VALUE=numeric(I),function(j) {
      rowSums(Y[,batch==j])/rowSums(Offset.new[,batch==j])
    })
    alphas <- log(mu)
    W <- matrix(0,nrow=I,ncol=J)
    for(j in 1:nbatch) {
      W[,batch==j] <- matrix(mu[,j], I, sum(batch==j))
    }
    W <- W*Offset.new
    W[W > max(Y)] <- max(Y)
    Z <- X+(Y-W)/W

    ## Compute log likelihood
    #LL[i] <- sum(dpois(Y,lambda=W,log=TRUE))
    LL[i] <- sum(dpois(Y,lambda=W,log=TRUE))
    cat("Iteration: ", i, ". LogLik=", LL[i], "\n")
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
    LRA <- irlba::irlba((W/w.max)*Z+(1-(W/w.max))*V,nv=M)
    X <- LRA$u %*% diag(LRA$d) %*% t(LRA$v)
  }

  out <- list()
  out$W <- W
  out$V <- LRA$v
  out$D <- LRA$d
  out$U <- LRA$u
  out$alpha <- alphas
  out$M <- M
  out$LL <- LL

  out <- process.results(out)
  return(out)
}

gbm.projection <- function(Y,M,subsample=2000,min.counts=5,ncores) {
  J <- ncol(Y)
  jxs <- sample(1:J,size=subsample,replace=FALSE)
  Y.sub <- Y[,jxs]
  Y.sub <- as.matrix(Y.sub)
  ixs <- which(rowSums(Y.sub) > 5)
  Y.sub <- Y.sub[ixs,]
  out <- gbm.sc(Y.sub,M=M)

  V.subsamp <- t(diag(out$D) %*% t(out$V))

  U <- out$U
  U <- as.data.frame(U)
  colnames(U) <- paste0("U",1:M)
  alpha <- out$alpha

  Y <- Y[ixs,]
  #Now get V
  V <- apply(Y,2,function(cell) {
    o <- log(sum(cell))+alpha
    fit <- glm(cell~0+offset(o)+.,
               data=U,
               family=poisson(link="log"))
    sumfit <- summary(fit)
    c(coefficients(fit),sumfit$coefficients[,2])
  })
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
  V <- t(diag(out$D)%*%t(out$V))
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

