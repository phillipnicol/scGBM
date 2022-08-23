

simData <- function(I,J,d,alpha.mean=0) {
  U <- rustiefel(m=I,R=d); V <- rustiefel(m=J,R=d)
  U <- U[,1:d]; V <- V[,1:d]
  D <- seq(2*(sqrt(I)+sqrt(J)),sqrt(I)+sqrt(J),length.out=d)
  for(m in 1:d) {
    if(U[1,m] < 0) {
      U[,m] <- -1*U[,m]
      V[,m] <- -1*V[,m]
    }
  }
  UDV <- U %*% diag(D) %*% t(V)
  alpha <- rnorm(n=I,mean=alpha.mean)
  Lambda <- matrix(alpha,nrow=I,ncol=J)+UDV
  Y <- matrix(rpois(n=I*J,lambda=as.vector(exp(Lambda))),nrow=I,ncol=J)

  val <- list()
  val$Y <- Y
  val$V <- t(diag(D) %*% t(V))
  val$U <- U
  val$D <- D
  return(val)
}










