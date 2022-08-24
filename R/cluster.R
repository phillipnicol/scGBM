



gbm.cluster <- function(out, se_V, max.cluster=100) {
  M <- length(out$D)
  V <- out$V
  J <- nrow(V)
  ll.null <- sum(dnorm(V/se_V,log=TRUE))

  bics <- rep(NA,max.cluster)
  bics[1] <- log(J)-2*ll.null

  for(K in 2:max.cluster) {
    print(K)
    fit <- kmeans(V,centers=K)

    Pi <- table(fit$cluster)/J
    mu <- fit$centers
    ll <- vapply(1:J,FUN.VALUE=numeric(1),function(j) {
      vec <- vapply(1:K,FUN.VALUE=numeric(1),function(k) {
        log(Pi[k])+sum(dnorm((V[j,]-mu[k,])/se_V[j,],log=TRUE))
      })
      c <- max(vec)
      c+log(sum(exp(vec-c)))
    })
    ll <- sum(ll)
    bics[K] <- (M*K+K)*log(J)-2*ll
    print(bics[K])

    if(bics[K-1] < bics[K]) {
      break
    }

    fit.prev <- fit

  }

  val <- list()
  val$bic <- bics
  val$cluster <- fit.prev$cluster
  val$fc.centers <- out$U %*% t(fit.prev$centers)
  A <- matrix(out$alpha,nrow=nrow(out$W),ncol=ncol(val$fc.centers))
  val$full.centers <- A + val$fc.centers
  return(val)
}

