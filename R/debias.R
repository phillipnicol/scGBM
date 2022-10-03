



debias.V <- function(out,Y) {
  I <- nrow(out$W); J <- ncol(out$W)
  M <- out$M
  Z <- out$U %*% t(out$V) + (Y-out$W)/out$W
  rWZ <- sqrt(out$W)*Z

  LRA <- svd(rWZ)

  sigma <- LRA$d
  w <- rep(0,length(sigma))
  for(m in 1:M) {
    if(sigma[m] > sqrt(I) + sqrt(J)) {
      w[m] <- 1-1/sigma[m]^2*(1+abs(J-I)+2*sum(sigma[m]^2/(sigma[m]^2-sigma[-m]^2)))
    }
  }

  V.scaled <- t(diag(1/out$D) %*% t(out$V))
  out$D <- out$D*w[1:M]
  out$D <- ifelse(out$D < 0, 0, out$D)
  out$V <- t(diag(out$D) %*% t(V.scaled))
  return(out)
}


robustify.V <- function(out) {
  V <- mclapply(1:J,function(j) {
    cell <- as.vector(Y[,j])
    o <- log(sum(cell))+alpha
    val <- matrix(rep(0,2*M),nrow=1)
    try({
      fit <- fastglm(x=U,y=cell,offset=o,
                     family=quasipoisson(),
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
