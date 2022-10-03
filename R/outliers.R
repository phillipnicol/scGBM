


find.outliers <- function(out,Y) {
  V <- data.frame(out$V)
  I <- nrow(Y); J <- ncol(Y)
  U <- out$U
  alpha <- out$alpha
  csy <- colSums(Y)
  C <- vapply(1:I,FUN.VALUE=numeric(J),function(i) {
    gene <- as.vector(Y[i,])
    o <- csy*exp(alpha[i])
    val <- rep(0,J)
    try({
      fit <- glm(gene~0+.,data=V,
                 family=poisson(),
                 offset=o)
      if(fit$converged) {
        val <- cooks.distance(fit)
      }
    })
    return(val)
  })


  V <- mclapply(1:J,function(j) {
    cell <- as.vector(Y[,j])
    o <- log(sum(cell))+alpha
    val <- matrix(rep(0,2*M),nrow=1)
    ixs <- which(C[j,]==1)
    try({
      fit <- fastglm(x=U[-ixs,],y=cell[-ixs],offset=o[-ixs],
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
}
