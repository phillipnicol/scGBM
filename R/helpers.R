

glmpcaObjective <- function(glmpca.fit,Y) {
  fit <- glmpca.fit
  I <- nrow(Y); J <- ncol(Y)

  U <- as.matrix(fit$loadings)
  V <- as.matrix(fit$factors)
  Alpha <- matrix(unlist(fit$coefX),nrow=I,ncol=J)
  O <- matrix(exp(fit$offsets), nrow=I, ncol=J, byrow=TRUE)

  Mu <- O*exp(Alpha+U%*%t(V))
  sum(Y*log(Mu)-Mu)
}

data.split <- function(Y,p) {
  I <- nrow(Y); J <- ncol(Y)
  Y1 <- matrix(rbinom(n=I*J, size=Y, prob=p),nrow=I,ncol=J)
  Y2 <- Y-Y1
  out <- list()
  out$Y1 <- Y1; out$Y2 <- Y2
  return(out)
}
