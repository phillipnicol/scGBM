
library(Seurat)
library(scGBM)

Y <- readRDS("../../data/blish_counts.RDS")  
Y <- as.matrix(Y)

set.seed(1)
subset <- seq(1000,43000,by=1000)
iters <- 25
Rmse <- matrix(0,nrow=length(subset),ncol=iters)

library(fastglm)
library(scGBM)
library(rstiefel)
library(doParallel)

out <- gbm.sc(Y,M=20)

true.v <- out$V


set.seed(1)
subset <- seq(1000,43000,by=1000)
iters <- 100
res <- array(dim=c(length(subset),iters,20))
for(i in 1:length(subset)) {
  cat("OVERALL ITERATION", i, " ", j, "\n")
  for(j in 1:iters) {
    out <- gbm.sc(Y,M=20,subset=subset[i],ncores=10)
    for(m in 1:20) {
    	  res[i,j,m] <- abs(cor(true.v[,m], out$V[,m]))
	  cat("RESULT:", res[i,j,m], "\n")
    }
  }
}


saveRDS(res, "../data/cor.RDS")