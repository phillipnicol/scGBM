



gbm.sc <- function(Y,M,max.iter=100,gamma.start=25,tol=10^{-4},
                   offset=NULL,subset=NULL,ncores=1) {
  if(!is.null(subset)) {
    if(ncores==1) {
      out <- gbm.projection(Y,M,subsample=subset)
    } else{
      out <- gbm.proj.parallel(Y,M,subsample=subset,ncores=ncores)
    }
    return(out)
  }

  if(is.null(offset)) {
    offset <- colSums(Y)
  }

  if(!is.matrix(Y)) {
    Y <- as.matrix(Y)
  }

  I <- nrow(Y); J <- ncol(Y)
  LL <- rep(0,max.iter)
  gamma <- gamma.start

  #Starting estimate for mu
  mu <- rowSums(Y)/sum(offset)
  alphas <- log(mu)
  Offset <- matrix(offset,nrow=I,ncol=J,byrow=TRUE)
  W <- matrix(mu,nrow=I,ncol=J)
  W <- W*Offset

  Z <- (Y-W)/sqrt(W)
  LRA <-  irlba::irlba(Z,nv=M,nu=M)
  X <- LRA$u %*% diag(LRA$d) %*% t(LRA$v)
  X <- sqrt(1/W)*X

  #Bound X to avoid starting too large
  X[X > 8] <- 8
  X[X < -8] <- -8

  last.change <- 1

  for(i in 1:max.iter) {
    #Reweight
    Offset.new <- Offset*exp(X)
    mu <- rowSums(Y)/rowSums(Offset.new)
    alphas <- log(mu)
    W <- matrix(mu,nrow=I,ncol=J)
    W <- W*Offset.new
    Z <- X+(Y-W)/W

    ## Compute log likelihood
    LL[i] <- sum(dpois(Y,lambda=W,log=TRUE))
    cat("Iteration: ", i, ". LogLik="
        , LL[i], ". Step size = ", gamma, ".\n")
    if(gamma < 0.01) {
      break
    }
    if(i > 1) {
      tau <- abs((LL[i]-LL[i-1])/LL[i-1])
      if(tau < tol) {
        break
      }

      if(LL[i-1] > LL[i]) {
        X <- X.prev
        gamma <- gamma/2
        last.change <- 0
        next
      } else {last.change <- last.change+1}

      if(last.change > 10) {
        gamma <- gamma*2
      }
    }

    ## Gradient Step
    X.prev <- X
    w.max <- max(W)
    LRA <- irlba::irlba(X-gamma*(W/w.max)*(X-Z),nv=M)
    X <- LRA$u %*% diag(LRA$d) %*% t(LRA$v)
  }

  out <- list()
  out$W <- W
  out$V <- LRA$v
  out$D <- LRA$d
  out$U <- LRA$u
  out$alpha <- alphas
  out$M <- M

  out <- process.results(out)
  return(out)
}

gbm.projection <- function(Y,M,subsample=2000,min.counts=5,ncores) {
  J <- ncol(Y)
  jxs <- sample(1:J,size=subsample,replace=FALSE)
  Y.sub <- Y[,jxs]
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
  jxs <- sample(1:J,size=subsample,replace=FALSE)
  Y.sub <- Y[,jxs]
  ixs <- which(rowSums(Y.sub) > 5)
  Y.sub <- Y.sub[ixs,]
  out <- gbm.sc(Y.sub,M=M,gamma.start = 50)

  U <- out$U
  U <- as.data.frame(U)
  colnames(U) <- paste0("U",1:M)
  alpha <- out$alpha

  Y <- Y[ixs,]
  #Now get V
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  V <- foreach(j=1:J, .combine=rbind) %dopar% {
    cell <- Y[,j]
    o <- log(sum(cell))+alpha
    fit <- glm(cell~0+offset(o)+.,
               data=U,
               family=poisson(link="log"))
    sumfit <- summary(fit)
    matrix(c(coefficients(fit),sumfit$coefficients[,2]),
              nrow=1)
  }

  out$se_V <- V[,(M+1):(2*M)]
  out$V <- V[,1:M]

  stopCluster(cl)
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



