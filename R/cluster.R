



gbm.cluster <- function(V, se_V, max.cluster=100) {
  M <- length(out$D)

  J <- nrow(V)
  ll.null <- ll
  ll.null <- sum(dnorm(V/se_V,log=TRUE))

  bics <- rep(NA,max.cluster)
  bics[1] <- log(J)-2*ll.null

  K <- max.cluster
  fit <- kmeans(V,centers=K)
  Pi <- table(fit$cluster)/J
  mu <- fit$centers
  ll <- 0
  ll <- vapply(1:J,FUN.VALUE=numeric(1),function(j) {
    vec <- vapply(1:K,FUN.VALUE=numeric(1),function(k) {
      log(Pi[k])+sum(dnorm((V[j,]-mu[k,])/se_V[j,],log=TRUE))
    })
    c <- max(vec)
    c+log(sum(exp(vec-c)))
  })
  ll <- sum(ll)
  bics[max.cluster] <- (M*K+K)*log(J)-2*ll

  K <- floor(max.cluster/2)
  lb <- 1; ub <- max.cluster
  K.prev <- 1
  while(TRUE) {
    print(K)
    fit <- kmeans(V,centers=K)

    if(K == max.cluster | K==1) {
      break
    } else if(abs(lb-ub) <= 1) {
      break
    }


    Pi <- table(fit$cluster)/J
    mu <- fit$centers
    ll <- 0
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

    if(bics[K] < bics[ub]) {
      ub <- K
      K.new <- floor(mean(c(ub,lb)))
    } else{
      lb <- K
      K.new <- floor(mean(c(ub,lb)))
    }
    K <- K.new
  }

  best <- which.min(bics)
  fit <- kmeans(V,centers=best)

  for(K in 2:max.cluster) {
    print(K)
    fit <- kmeans(V,centers=K)
    Pi <- table(fit$cluster)/J
    mu <- fit$centers
    ll <- 0
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

    if(bics[K] > bics[K-1]) {
      break
    }
    cluster <- fit$cluster
  }

  val <- list()
  val$bic <- bics
  val$cluster <- cluster
  return(val)
}


