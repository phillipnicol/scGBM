
library(DuoClustering2018)
library(scGBM)

sce <- sce_full_Zhengmix8eq() 

Y <- counts(sce)
Y <- as.matrix(Y)
Y <- Y[rowSums(Y) >= 10,]
set.seed(1)
J <- ncol(Y)
#perm <- sample(1:J,size=J,replace=FALSE)
#Rmse <- matrix(0,nrow=length(subset),ncol=iters)

library(fastglm)
library(scGBM)
library(rstiefel)
library(doParallel)

out <- gbm.sc(Y,M=20)

true.v <- out$V |> as.matrix()

Pv.true <- out$V %*% t(out$V)

Pv.true.top5 <- out$V[,1:5] %*% t(out$V[,1:5])

print("STARTING")
set.seed(1)
#subset <- seq(100,1000,by=100)
#subset <- c(subset, 3994)
subset <- seq(100, 3994, length.out=15)
iters <- 10
res <- array(dim=c(length(subset),iters,20))
res.dist <- matrix(0, nrow=length(subset), ncol=iters)
res.dist5 <- matrix(0, nrow=length(subset), ncol=iters)

for(i in 1:length(subset)) {
  for(j in 1:iters) {
    out <- gbm.sc(Y,M=20,subset=subset[i],ncores=10)
    out$V <- as.matrix(out$scores)
    cat("OVERALL ITERATION", i, " ", j, "\n")
    for(m in 1:20) {
    	  res[i,j,m] <- abs(cor(true.v[,m],out$V[,m]))
          cat("RESULT:", res[i,j,m], "\n")
	  
    }

    Pv.est <- out$V %*% solve(t(out$V) %*% out$V) %*% t(out$V)

    V.hat5 <- out$V[,1:5]
    Pv.top5.est <- V.hat5 %*% solve(t(V.hat5) %*% V.hat5) %*% t(V.hat5)

    res.dist[i,j] <- sqrt(sum((Pv.true - Pv.est)^2))
    res.dist5[i,j] <- sqrt(sum((Pv.true.top5 - Pv.top5.est)^2))
  }
}


saveRDS(res, "../data/corZheng.RDS")
saveRDS(res.dist, "../data/proj_dist_Zheng.RDS")
saveRDS(res.dist, "../data/proj_dist_Zheng_top5.RDS")