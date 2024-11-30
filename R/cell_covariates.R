#' @export
#'
#' @title scGBM: Model-based dimensionality reduction for single-cell
#' RNA-seq with generalized bilinear models
#'
#' @description Fit a Poisson bilinear model to the single-cell
#' count matrix.
#'
#' @param Y A matrix of UMI counts with genes on rows and cells on columns.
#' @param M The number of latent factors to estimate
#' @param max.iter The maximum number of iterations
#' @param tol The tolerance for convergence (relative difference in objective function)
#' @param subset If NULL, use the entire data to fit the model. If an integer, take a
#' random subsample of cells to fit the loadings. Then project the remaining
#' cells onto those loadings. If a non-random sample is desired, subset can
#' be an integer vector corresponding to the columns that should be used to
#' fit the loadings.
#' @param ncores If subset is not NULL, then ncores specifies the number of cores to use
#' during the projection step.
#' @param infer.beta If FALSE, beta is set to the log of the number of counts per cell.
#' @param return.W If TRUE, the matrix of weights (which is equal to the estimated mean) is returned. This should
#' be FALSE for large datasets since it requires a lot of memory.
#' @param batch An optional factor containing the assignment of cells to known batches.
#' @param time.by.iter If TRUE, the elapsed time (in seconds) is given at each iteration of the algorithm.
#' @param min.iter The minimum number of iterations.
#'
#' @return A list with components
#' \itemize{
#' \item \code{scores} - A matrix containing the factor scores (weight of each cell for each factor).
#' \item \code{loadings} - A matrix (genes x factors) containing the factor loadings (weight of each cell for each factor).
#' \item \code{V} - A matrix containing the unscaled factor scores.
#' \item \code{U} - Factor loadings. This matrix is identical to \code{loadings} and is included
#' for reverse compatibility.
#' \item \code{D} - A vector containing the singular values (scaling factors)
#' \item \code{alpha} - A matrix containing gene-specific intercepts. The number of columns
#' is set to be equal to the number of batches (by default there are no batches so this is 1).
#' \item \code{beta} - A vector of cell-specific intercepts.
#' \item \code{I} - The number of genes.
#' \item \code{J} - The number of cells.
#' \item \code{W} - The estimated mean (and by properties of the Poisson distribution, also the variance)
#' for each entry of the count matrix.
#' \item \code{obj} - The value of the objective function for each iteration.
#' }
#'
#' @details scGBM fits the following model to the \eqn{I \times J} count matrix \eqn{Y}:
#' \deqn{Y \sim \text{Poisson}(\mu)}
#' \deqn{\log \mu = \alpha 1_J^T + 1_I \beta^T + U \Sigma V^T}
#' Here \eqn{\alpha} are gene-specific intercepts, \eqn{\beta_j} are cell-specific intercepts, \eqn{U} are the \eqn{I \times M} factor loadings
#' (linear combinations of genes that define factors), and \eqn{V \Sigma} are the \eqn{J \times M} scores for each cell.
#'
#' @author Phillip B. Nicol <philnicol740@gmail.com>
gbm.cell.covariate <- function(Y,
                               M,
                               max.iter=100,
                               tol=10^{-4},
                               subset=NULL,
                               ncores=1,
                               infer.beta=FALSE,
                               return.W = TRUE,
                               X.cov = matrix(1,nrow=ncol(Y), ncol=1),
                               time.by.iter = FALSE,
                               min.iter=30) {

  #Check validity of input
  gbm.sc.check.valid.input(as.list(environment()))

  if(!is.null(subset)) {
    out <- gbm.proj.parallel(Y,M,subsample=subset,ncores=ncores,tol=tol,
                             max.iter=max.iter)
    ##Message for users of new version about scores
    message("For users of newer versions (1.0.1+): the `scores` matrix now contains factor scores.")
    return(out)
  }


  I <- nrow(Y); J <- ncol(Y)
  LL <- rep(0,max.iter)
  loglik <- c()
  if(time.by.iter) {
    time <- c()
    start.time <- Sys.time()
  }

  lr <- 1 #Default learning rate is set to 1


  #Precompute relevant quantities
  max.Y <- max(Y)
  nz <- which(Y != 0)
  log.csY <- log(colSums(Y))

  #Starting estimate for alpha and W
  betas <- log.csY

  W <- matrix(0, nrow=I, ncol=J)

  A <- matrix(nrow=I,ncol=ncol(X.cov))
  for(ix in 1:I) {
    A[ix,] <- coefficients(glm(Y[ix,] ~ . + 0, data=as.data.frame(X.cov),
                              offset = betas,
                              family=poisson()))
  }



  #betas <- betas + mean(alphas)
  #alphas <- alphas - mean(alphas) #Ensure alphas sum to 0
  W <- exp(sweep(A %*% t(X.cov), 2, betas, "+"))

  #Starting estimate of X
  Z <- (Y-W)/sqrt(W)
  LRA <-  irlba::irlba(Z,nv=M,nu=M)
  X <- LRA$u %*%(LRA$d*t(LRA$v))
  X <- sqrt(1/W)*X

  X[X > 8] <- 8 #Clipping
  X[X < -8] <- -8

  #For acceleration, save previous X
  Xt <- matrix(0,nrow=I,ncol=J)

  for(i in 1:max.iter) {
    #Reweight
    A <- matrix(nrow=I,ncol=ncol(X.cov))
    for(ix in 1:I) {
      A[ix,] <- coefficients(glm(Y[ix,] ~ . + 0, data=as.data.frame(X.cov),
                                offset = betas + X[ix,],
                                family=poisson()))
    }

    #betas <- betas + mean(alphas)
    #alphas <- alphas - mean(alphas)
    if(infer.beta) {
      betas <- log.csY-log(colSums(exp(alphas[,batch]+X)))
    }

    W <- exp(sweep(A %*% t(X.cov) + X, 2, betas, "+"))

    ## Timing
    if(time.by.iter) {
      time <- c(time,difftime(Sys.time(),start.time,units="sec"))
      start.time <- Sys.time()
    }

    #Prevent W from being too large (stability)
    W[W > max.Y] <- max.Y

    #Compute working variables
    #Z <- X+(Y-W)/W

    ## Compute log likelihood (no normalizing constant)
    LL[i] <- sum(Y[nz]*log(W[nz]))-sum(W)
    if(is.na(LL[i]) | is.infinite(LL[i])) {
      X <- Xt
      lr <- lr/2
      i <- i - 1
      next
    }
    if(i >= 3) {
      tau <- abs((LL[i]-LL[i-2])/LL[i])
      if(tau < tol & lr <= 1.06 & i >= min.iter) {
        break
      }

      if(LL[i] < LL[i-1]) {
        lr <- max(lr/2, 1)
        X <- Xt
        next
      } else {
        lr <- lr*(1.05)
        #lr <- lr
      }
    }

    loglik <- c(loglik,LL[i])
    cat("Iteration: ", i, ". Objective=", LL[i], "\n")

    ### Projected gradient descent step
    pgd <- pgd_irlba(X, Xt, i, lr, W, Y, M)
    X <- pgd$X
    Xt <- pgd$Xt

    if(i == max.iter) {
      warning("Maximum number of iterations reached (increase max.iter).
              Possible non-convergence.")
    }
  }

  out <- list()
  if(return.W) {
    out$W <- exp(sweep(A %*% t(X.cov) + X, 2, betas, "+"))
  }
  LRA <- pgd$LRA
  out$V <- LRA$v; rownames(out$V) <- colnames(Y); colnames(out$V) <- 1:M
  out$D <- LRA$d; names(out$D) <- 1:M
  out$U <- LRA$u; rownames(out$U) <- rownames(Y); colnames(out$U) <- 1:M
  out$loadings <- out$U
  out$A <- A
  out$beta <- betas; names(out$beta) <- colnames(Y)
  out$M <- M
  out$I <- nrow(out$W); out$J <- ncol(out$W)
  out$obj <- loglik
  if(time.by.iter) {
    out$time <- cumsum(time)
  }

  out <- process.results(out)

  ##Message for users of new version about scores
  message("For users of newer versions (1.0.1+): the `scores` matrix now contains factor scores, the `V` matrix is UNSCALED scores.")
  return(out)
}
