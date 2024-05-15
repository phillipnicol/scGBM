library(rsteifel)
library(glmpca)
library(scGBM)


I <- 10^3
J <- 10^4
M <- 10

reps <- 100

#Fisher, Avagrad, SGD, scGBM, scGBM-proj
nmethod <- 5

res <- array(0,dim=c(reps,nmethod,M)) #For factor cor
res2 <- matrix(0,nrow=reps, ncol=nmethod) #For X mse

find.best.cor <- function(V.true, factor.test) {
  val <- rep(0, ncol(V.true))
  for(m in 1:ncol(V.true)) {
    val[m] <- abs(cor(V.true[,m], factor.test))
  }
  return(max(val))
}

set.seed(1)
for(i in 1:reps) {
  sim <- simData(I,J,d=M)
  X.true <- sim$U %*% t(sim$V)

  #gbm-sc
  out <- gbm.sc(sim$Y,M=10,order.by.deviance = FALSE)
  for(m in 1:M) {
    res[i,1,m] <- find.best.cor(sim$V,out$scores[,m])
  }
  X.pred <- out$U %*% t(out$scores)
  res2[i,1] <- sqrt(mean((X.true - X.pred)^2))

  out <- gbm.sc(sim$Y,M=10,order.by.deviance = FALSE,
                subset=10^3, ncores=8)
  for(m in 1:M) {
    res[i,2,m] <- find.best.cor(sim$V,out$scores[,m])
  }
  X.pred <- out$U %*% t(out$scores)
  res2[i,2] <- sqrt(mean((X.true - X.pred)^2))

  fit <- glmpca(Y=sim$Y,L=10)
  fit$res$factors <- as.matrix(fit$res$factors)
  fit$res$loadings <- as.matrix(fit$res$loadings)
  for(m in 1:M) {
    res[i,3,m] <- find.best.cor(sim$V,fit$res$factors[,m])
  }
  X.pred <- fit$res$loadings %*% t(fit$res$factors)
  res2[i,3] <- sqrt(mean((X.true - X.pred)^2))

  fit <- glmpca(Y=sim$Y,L=10,minibatch="stochastic",ctl=list(batch_size=200))
  fit$res$factors <- as.matrix(fit$res$factors)
  fit$res$loadings <- as.matrix(fit$res$loadings)
  for(m in 1:M) {
    res[i,4,m] <- find.best.cor(sim$V,fit$res$factors[,m])
  }
  X.pred <- fit$res$loadings %*% t(fit$res$factors)
  res2[i,4] <- sqrt(mean((X.true - X.pred)^2))

  fit <- glmpca(Y=sim$Y, L=10,optimizer = "fisher")
  fit$res$factors <- as.matrix(fit$res$factors)
  fit$res$loadings <- as.matrix(fit$res$loadings)
  for(m in 1:M) {
    res[i,5,m] <- find.best.cor(sim$V,fit$res$factors[,m])
  }
  X.pred <- fit$res$loadings %*% t(fit$res$factors)
  res2[i,5] <- sqrt(mean((X.true - X.pred)^2))
}

saveRDS(res, "../data/factor_correlation.RDS")
saveRDS(res2, "../data/mse.RDS")

