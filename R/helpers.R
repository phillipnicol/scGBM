

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
