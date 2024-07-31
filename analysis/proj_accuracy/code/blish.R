

#library(Seurat)
#Sco <- readRDS("../../data/blish.RDS")

Y <- readRDS("../../data/blish_counts.RDS")
Y <- as.matrix(Y)

set.seed(1)
subset <- seq(1000,10000,by=1000)
iters <- 25
Rmse <- matrix(0,nrow=length(subset),ncol=iters)

library(fastglm)
library(scGBM)
library(rstiefel)
library(doParallel)

out <- gbm.sc(Y,M=20,max.iter=250,tol=0)

true.v <- out$V
Pv <- out$V %*% t(out$V)

set.seed(1)
subset <- seq(1000,43000,by=1000)
iters <- 100
res <- matrix(0, nrow=length(subset), ncol=iters)
for(i in 1:length(subset)) {
  for(j in 1:iters) {
    out <- gbm.sc(Y,M=20,subset=subset[i],ncores=10)
    Pv.hat <- out$scores %*% solve(t(out$scores) %*% out$scores + diag(0.001)) %*% t(out$scores)
    res[i,j] <- sqrt(mean((Pv - Pv.hat)^2))
    print(res[i,j])
    saveRDS(res, "blish_proj_accuracy.RDS")
  }
}
