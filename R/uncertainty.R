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
get.se <- function(out, EPS=0) {
  M <- length(out$D)
  U <- out$U
  V <- out$scores
  W <- out$W
  J <- nrow(V)
  tU <- t(out$U)

  se_V <- apply(W,2,function(w) {
    Mat <- tU %*% (w*out$U)
    sqrt(diag(solve(Mat + diag(EPS,nrow=M))))
  })
  out$se_scores <- t(se_V)

  tV <- t(V)
  se_U <- apply(W,1,function(w) {
    Mat <- (tV %*% (w*V))
    sqrt(diag(solve(Mat + diag(EPS,nrow=M))))
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

  U.denoise <- out$loadings
  for(m in 1:M) {
    U.denoise[,m] <- ashr::get_pm(ashr::ash(U.denoise[,m],
                                            out$se_loadings[,m]))
  }

  out$loadings.denoised <- U.denoise
  return(out)
}







get.se2 <- function(out, EPS=0) {
  M <- length(out$D)
  U <- out$U
  V <- out$scores
  W <- out$W
  J <- nrow(V)
  tU <- t(out$U)

  US <- out$U %*% diag(out$D)
  tUS <- t(US)

  n <- 100
  U.post <- array(U,dim=c(I,M,n))

  V.post <- array(0, dim=c(J,M,n))
  for(k in 1:n) {
    V.post[,,k] <- t(apply(W,2,function(w) {
      Mat <- (t(U.post[,,k]) %*% (w*U.post[,,k]))
      mvtnorm::rmvnorm(n=1,sigma=solve(Mat))
    }))
    V.post[,,k] <- V.post[,,k] + V
    V.post[,,k] <- normalize.cols(V.post[,,k])
    V.post[,,k] <- V.post[,,k] %*% diag(out$D)
  }

  sdV <- aaply(V.post,1:2,sd)
  print(mean(sdV))

  U.post <- array(0, dim=c(I,M,n))
  for(k in 1:n) {
    U.post[,,k] <- t(apply(W,1,function(w) {
      Mat <- (t(V.post[,,k]) %*% (w*V.post[,,k]))
      mvtnorm::rmvnorm(n=1,sigma=solve(Mat))
    }))
    U.post[,,k] <- U.post[,,k] + U
    U.post[,,k] <- normalize.cols(U.post[,,k])
  }
  sdU <- aaply(U.post,1:2,sd)
  print(mean(sdU))

  tV <- t(V)
  se_U <- apply(W,1,function(w) {
    Mat <- (tV %*% (w*V))
    sqrt(diag(solve(Mat + diag(EPS,nrow=M))))
  })
  out$se_loadings <- t(se_U)
  return(out)
}



