

simData <- function(I,J,d,alpha.mean=0,beta.mean=0) {
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
  beta <- rnorm(n=J,mean=beta.mean)
  Lambda <- matrix(alpha,nrow=I,ncol=J)+matrix(beta,nrow=I,ncol=J,byrow=TRUE)+UDV
  Y <- matrix(rpois(n=I*J,lambda=as.vector(exp(Lambda))),nrow=I,ncol=J)

  val <- list()
  val$Y <- Y
  val$V <- t(diag(D) %*% t(V))
  val$U <- U
  val$D <- D
  val$alpha <- alpha
  val$beta <- beta
  return(val)
}


simNull <- function(I,J) {
  mus <- runif(n=I,min=1,max=10)
  Mu <- matrix(mus,nrow=I,ncol=J)

  Y <- matrix(rpois(n=I*J,lambda=as.vector(Mu)),nrow=I,ncol=J)
  return(Y)
}


simNullNB <- function(I,J,theta=1) {
  mus <- runif(n=I,min=1,max=10)
  Mu <- matrix(mus,nrow=I,ncol=J)
  Y <- matrix(MASS::rnegbin(n=I*J,mu=as.vector(Mu),theta=theta),nrow=I,ncol=J)
  return(Y)
}



