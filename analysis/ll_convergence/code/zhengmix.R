
library(scGBM)
library(rstiefel)

#Compute log-likelihood
gbm.ll <- function(W,Y) {
  sum(dpois(Y,lambda=W,log=TRUE))
}
glmpca.ll <- function(U,V,alpha,beta,Y) {
  I <- nrow(Y); J <- ncol(Y)

  Alpha <- matrix(alpha,nrow=I,ncol=J)
  O <- matrix(exp(beta), nrow=I, ncol=J, byrow=TRUE)

  Mu <- O*exp(Alpha+U%*%t(V))
  sum(dpois(Y,Mu,log=TRUE))
}


print(R.version)


set.seed(1)
library(DuoClustering2018)
library(Seurat)
library(scGBM)
library(fastglm)
library(bigmemory)



sce <- sce_full_Zhengmix8eq()

Sco <- as.Seurat(sce)
Sco <- NormalizeData(Sco)
Sco <- FindVariableFeatures(Sco,nfeatures=2000)
Y <- Sco@assays$originalexp@counts[Sco@assays$originalexp@var.features,]
Y <- as.matrix(Y)

Y <- Y[rowSums(Y) >= 50,] 
ds <- scGBM:::data.split(Y,p=0.5)
Y1 <- ds$Y1; Y2 <- ds$Y2
I <- nrow(Y); J <- ncol(Y)

max.iter <- 250
out <- gbm.sc(Y1,oos.Y=Y2,M=20,max.iter=max.iter,tol=10^{-5},infer.beta=TRUE,time.by.iter = TRUE)
print(out$ll.oos)

time.1 <- out$time
ll.1 <- out$ll.oos

## GBM SC PAR
## GBM SC PAR
library(parallel)
library(fastglm)

proj.time <- c()
proj.ll <- c()
subset <- sample(1:ncol(Y),size=400,replace=FALSE)
for(k in seq(25,250,by=25)) {
  print(k)
  start <- Sys.time()
  out <- gbm.sc(Y1,M=20,max.iter=k,subset=subset,ncores=12,tol=10^{-4})
  end <- Sys.time()
  proj.time <- c(proj.time,difftime(end,start,units="secs")[[1]])

  Alpha <- matrix(out$alpha,nrow=I,ncol=J)
  Beta <- matrix(out$beta,nrow=I,ncol=J,byrow=TRUE)
  W <- exp(Alpha+Beta+out$U %*% t(out$scores))
  proj.ll <- c(proj.ll,sum(Y2*log(W) - W))
  print(sum(Y2*log(W) - W))
}
ll.2 <- proj.ll
time.2 <- proj.time





## GLM-PCA Avagrad
library(glmpca)
max.iter <- 250

time.3 <- c(1:max.iter)
ll.3 <- c(1:max.iter)

fit <- glmpca(Y1,L=20,Y.oos=Y2,ctl=list(maxIter=max.iter))
time.3 <- fit$mylist$time[-1]
ll.3 <- fit$mylist$LL

max.iter <- 100
fit <- glmpca(Y1,L=20,Y.oos=Y2,optimizer="fisher",ctl=list(verbose=TRUE,maxIter=max.iter))


time.4 <- fit$mylist$time[-1]
ll.4 <- fit$mylist$LL

set.seed(1)
max.iter <- 500
fit <- glmpca(Y1,Y.oos=Y2,L=20,minibatch="stochastic",ctl=list(verbose=TRUE,max.iter=max.iter,batch_size=400))


time.5 <- fit$mylist$time[-1]
ll.5 <- fit$mylist$LL


nalgos <- 3
Time <- list(time.1,time.2,time.3,time.4, time.5)
LL <- list(ll.1,ll.2,ll.3,ll.4, ll.5)


save(Time,file="../data/zhengmix_time.RData")
save(LL,file="../data/zhengmix_LL.RData")