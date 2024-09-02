
library(DuoClustering2018)
library(scGBM)

sce <- sce_full_Zhengmix8eq() 

Y <- counts(sce)
Y <- as.matrix(Y)
set.seed(1)
J <- ncol(Y)
#perm <- sample(1:J,size=J,replace=FALSE)
#Rmse <- matrix(0,nrow=length(subset),ncol=iters)

library(fastglm)
library(scGBM)
library(rstiefel)
library(doParallel)

out <- gbm.sc(Y,M=20)

true.v <- out$V


set.seed(1)
subset <- seq(100,3900,by=100)
iters <- 10
res <- array(dim=c(length(subset),iters,20))
for(i in 1:length(subset)) {
  for(j in 1:iters) {
    out <- gbm.sc(Y,M=20,subset=subset[i],ncores=10)
    cat("OVERALL ITERATION", i, " ", j, "\n")
    for(m in 1:20) {
    	  res[i,j,m] <- abs(cor(true.v[,m],out$V[,m]))
          cat("RESULT:", res[i,j,m], "\n")
	  
    }
  }
}


saveRDS(res, "../data/corZheng.RDS")