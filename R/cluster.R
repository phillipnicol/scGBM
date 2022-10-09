


gbm.cluster2 <- function(out,se_V,max.cluster=20,keep) {
  M <- keep
  V <- out$V[,1:keep]
  J <- nrow(V)
  se_V <- se_V[,1:keep]

  fit <- Mclust(V,G=1,model="VVI")

  shrinkage <- colMeans(se_V^2)
  Sigma <- fit$parameters$variance$Sigma - diag(shrinkage)
  Sigma <- ifelse(Sigma < 0, 0, Sigma)
  ll <- vapply(1:J,FUN.VALUE=numeric(1),function(j) {
    sum(dnorm(V[j,],mean=t(fit$parameters$mean),
              sd=sqrt(diag(Sigma)+se_V[j,]^2),log=TRUE))
  })

  bics <- rep(NA,max.cluster)
  bics[1] <- 2*M*log(J)-2*sum(ll)
  fit.prev <- fit
  for(K in 2:max.cluster) {
    print(K)
    fit <- Mclust(V,G=K,model="VVI")
    mu <- t(fit$parameters$mean)
    Pi <- fit$parameters$pro

    Sigma <- list()
    w <- estep(data=V,modelName="VVI",parameters=fit$parameters)$z
    for(k in 1:K) {
      shrinkage <- apply(se_V^2,2,function(x) weighted.mean(x,w=w[,k]))
      Sigma[[k]] <- fit$parameters$variance$sigma[,,k]-diag(shrinkage)
      Sigma[[k]] <- ifelse(Sigma[[k]] < 0, 0, Sigma[[k]])
    }

    ll <- vapply(1:J,FUN.VALUE=numeric(1),function(j) {
      vec <- vapply(1:K,FUN.VALUE=numeric(1),function(k) {
        log(Pi[k])+sum(dnorm(V[j,],mean=mu[k,],
                             sd=sqrt(diag(Sigma[[k]])+se_V[j,]^2),
                             log=TRUE))
      })
      c <- max(vec)
      c+log(sum(exp(vec-c)))
    })
    ll <- sum(ll)
    bics[K] <- (2*M*K+K)*log(J)-2*ll
    print(bics[K])

    #if(bics[K-1] < bics[K]) {
    #  break
    #}

  }
  val <- list()
  #Apply evanno method
  if(bics[1] < bics[2]) {
    #K=1
    val$cluster <- rep(1,J)
  } else{
    #K=2
    cl <- which.min(bics)
    fit <- Mclust(V,G=cl,model="VVI")
    val$cluster <- fit$classification
  }

  #Shrinkage
  #mu <- t(fit.prev$parameters$mean)
  val$bic <- bics
  #val$fc.centers <- out$U[,1:keep] %*% mu
  #A <- matrix(out$alpha,nrow=nrow(out$U),ncol=ncol(val$fc.centers))
  #val$full.centers <- A + val$fc.centers
  return(val)
}


gbm.cluster3 <- function(out,se_V,max.cluster=20,keep) {
  wi <- 30
  M <- keep
  V <- out$V[,1:keep]
  J <- nrow(V)
  se_V <- se_V[,1:keep]

  fit <- Mclust(V,G=1,model="VVI")

  shrinkage <- colMeans(se_V^2)
  Sigma <- fit$parameters$variance$Sigma - diag(shrinkage)
  Sigma <- ifelse(Sigma < wi, wi, Sigma)
  ll <- vapply(1:J,FUN.VALUE=numeric(1),function(j) {
    sum(dnorm(V[j,],mean=t(fit$parameters$mean),
              sd=sqrt(diag(Sigma)+se_V[j,]^2),log=TRUE))
  })

  bics <- rep(NA,max.cluster)
  bics[1] <- 2*M*log(J)-2*sum(ll)
  fit.prev <- fit
  for(K in 2:max.cluster) {
    print(K)
    fit <- Mclust(V,G=K,model="VVI")
    mu <- t(fit$parameters$mean)
    Pi <- fit$parameters$pro

    Sigma <- list()
    w <- estep(data=V,modelName="VVI",parameters=fit$parameters)$z
    for(k in 1:K) {
      shrinkage <- apply(se_V^2,2,function(x) weighted.mean(x,w=w[,k]))
      Sigma[[k]] <- fit$parameters$variance$sigma[,,k]-diag(shrinkage)
      Sigma[[k]] <- ifelse(Sigma[[k]] < wi, wi, Sigma[[k]])
    }

    ll <- vapply(1:J,FUN.VALUE=numeric(1),function(j) {
      vec <- vapply(1:K,FUN.VALUE=numeric(1),function(k) {
        log(Pi[k])+sum(dnorm(V[j,],mean=mu[k,],
                             sd=sqrt(diag(Sigma[[k]])+se_V[j,]^2),
                             log=TRUE))
      })
      c <- max(vec)
      c+log(sum(exp(vec-c)))
    })
    ll <- sum(ll)
    bics[K] <- (2*M*K+K)*log(J)-2*ll
    print(bics[K])

    #if(bics[K-1] < bics[K]) {
    #  break
    #}

  }
  val <- list()
  #Apply evanno method
  if(bics[1] < bics[2]) {
    #K=1
    val$cluster <- rep(1,J)
  } else{
    #K=2
    cl <- which.min(bics)
    fit <- Mclust(V,G=cl,model="VVI")
    val$cluster <- fit$classification
  }

  #Shrinkage
  mu <- t(fit.prev$parameters$mean)
  val$bic <- bics
  #val$fc.centers <- out$U[,1:keep] %*% mu
  #A <- matrix(out$alpha,nrow=nrow(out$U),ncol=ncol(val$fc.centers))
  #val$full.centers <- A + val$fc.centers
  return(val)
}

gbm.cluster <- function(out, se_V, max.cluster=100) {
  M <- length(out$D)
  V <- out$V
  J <- nrow(V)
  ll.null <- sum(dt(V/se_V,log=TRUE,df=1))

  bics <- rep(NA,max.cluster)
  bics[1] <- log(J)-2*ll.null

  for(K in 2:max.cluster) {
    print(K)
    fit <- kmeans(x=V,centers=K,nstart=25)

    Pi <- table(fit$cluster)/J
    mu <- fit$centers
    ll <- vapply(1:J,FUN.VALUE=numeric(1),function(j) {
      vec <- vapply(1:K,FUN.VALUE=numeric(1),function(k) {
        log(Pi[k])+sum(dt((V[j,]-mu[k,])/se_V[j,],df=1,log=TRUE))
      })
      c <- max(vec)
      c+log(sum(exp(vec-c)))
    })
    ll <- sum(ll)
    bics[K] <- (M*K+K)*log(J)-2*ll
    print(bics[K])

    #if(bics[K-1] < bics[K]) {
    #  break
    #}

    fit.prev <- fit

  }

  val <- list()
  val$bic <- bics
  val$cluster <- fit.prev$cluster
  val$fc.centers <- out$U %*% t(fit.prev$centers)
  A <- matrix(out$alpha,nrow=nrow(out$U),ncol=ncol(val$fc.centers))
  val$full.centers <- A + val$fc.centers
  return(val)
}


sigma.shrink <- function(mclust.fit, K) {
  w <- seq(0,1,by=0.01)
  if(K == 1) {
    Sigma <- fit$parameters$variance$Sigma
    mu <- t(fit$parameters$mean)
    ll <- vapply(1:J,FUN.VALUE=numeric(1),function(j) {
      sum(dnorm(V[j,],mean=mu,sd=sqrt(se_V[j,]^2+w*diag(Sigma)),log=TRUE))
    })
  }
}

model.lik <- function(V,se_V,Sigma,K) {
  ll <- vapply(1:J,FUN.VALUE=numeric(1),function(j) {
    vec <- vapply(1:K,FUN.VALUE=numeric(1),function(k) {
      log(Pi[k])+sum(dnorm())
    })
    c <- max(vec)
    c+log(sum(exp(vec-c)))
  })
}
