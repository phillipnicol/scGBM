#' @export
#'
#' @title Uncertainty quantification
#'
#' @description Perform uncertainty quantification in the low-dimensional representation
#' estimated by scGBM.
#'
#' @param out a list that is the return value of \code{gbm.sc}
#'
#' @return Appends a two matrices \code{se_U} and \code{se_V} to
#' the list out containing the standard errors.
#'
#' @author Phillip B. Nicol <philnicol740@gmail.com>
get.se <- function(out) {
  M <- length(out$D)
  U <- out$U
  V <- out$V
  W <- out$W
  J <- nrow(V)
  tU <- t(out$U)
  EPS <- 0.001

  se_V <- apply(W,2,function(w) {
    Mat <- tU %*% (w*out$U)
    sqrt(diag(solve(Mat)))
  })
  out$se_scores <- t(se_V)

  tV <- t(V)
  se_U <- apply(W,1,function(w) {
    Mat <- (tV %*% (w*V))
    sqrt(diag(solve(Mat)))
  })
  out$se_loadings <- t(se_U)
  return(out)
}

#' @export
#'
#' @title Denoise loadings
#'
#' @description Use the adaptive shrinkage method implement in \code{ashr}
#' to stabilize estimates of U.
#'
#' @author Phillip B. Nicol <philnicol740@gmail.com>
denoise.U <- function(out) {
  M <- length(out$D)

  U.denoise <- out$U
  for(m in 1:M) {
    U.denoise[,m] <- ashr::get_pm(ashr::ash(U.denoise[,m],out$se_U[,m]))
  }

  out$loadings.denoised <- U.denoise
  return(out)
}
