#' @export
get.se <- function(out) {
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
  out$se_V <- t(se_V)

  tV <- t(V)
  se_U <- apply(W,1,function(w) {
    Mat <- (tV %*% (w*V))
    sqrt(diag(solve(Mat)))
  })
  out$se_U <- t(se_U)
  return(out)
}

denoise.U <- function(out) {
  M <- length(out$D)

  U.denoise <- out$U
  for(m in 1:M) {
    U.denoise[,m] <- ashr::get_pm(ashr::ash(U.denoise[,m],out$se_U[,m]))
  }

  out$U.denoise <- U.denoise
  return(out)
}


## Quantify uncertainty

## Quantify uncertainty

uncertainty <- function(out) {
  ## Get corU
  M <- length(out$D)
  U <- out$U
  V <- out$V
  W <- out$W
  J <- nrow(V)
  I <- nrow(U)

  tV <- t(V)
  cov.U <- t(array(apply(W,1,function(w) {
    Mat <- tV %*% (w*V)
    diag(solve(Mat))
  }),dim=c(M,I)))
  cov.U.init <- cov.U

  tU <- t(U)
  cov.V <- t(array(apply(W,2,function(w) {
    Mat <- tU %*% (w*U)
    diag(solve(Mat))
  }),dim=c(M,J)))
  cov.V.init <- cov.V

  #se.V <- sqrt(cov.V[1,1,])

  X <- U %*% t(V)
  Z <- X+(Y-W)/W
  W.scale <- W/max(W)
  Q <- W.scale*Z+(1-W.scale)*X
  Q2 <- Q^2

  iter.max <- 50
  for(r in 1:iter.max) {
    varVinvVtV <- matrix(0,nrow=J,ncol=M)
    for(m in 1:M) {
      sm <- sum(V[,m]^2)
      nabla.f <- -2*outer(V[,m],V[,m])/sm^2
      diag(nabla.f) <- diag(nabla.f) + 1/sm
      varVinvVtV[,m] <- diag(nabla.f %*% t(nabla.f * cov.V[,m]))
    }
    cov.U <- cov.U.init + Q2 %*% varVinvVtV

    cov.V <- cov.V.init + t(Q2) %*% cov.U
    print(sqrt(cov.V[1,1]))
  }

}

