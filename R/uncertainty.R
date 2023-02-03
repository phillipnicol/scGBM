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
