devtools::install_github("phillipnicol/scGBM", ref="dev2")


library(scGBM)
library(rstiefel)
library(SingleCellExperiment)
library(Seurat)

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
library(scGBM)
library(fastglm)
library(bigmemory)



sce <- sce_full_Zhengmix8eq()

#Sco <- as.Seurat(sce)
#Sco <- NormalizeData(Sco)
#Sco <- FindVariableFeatures(Sco,nfeatures=2000)
#Y <- Sco@assays$originalexp@counts
Y <- counts(sce)
Y <- as.matrix(Y)

Y <- Y[rowSums(Y) >= 50,]
ds <- scGBM:::data.split(Y,p=0.5)
Y1 <- ds$Y1; Y2 <- ds$Y2
I <- nrow(Y); J <- ncol(Y)
colnames(Y1) <- colnames(Y); rownames(Y1) <- rownames(Y)
print(dim(Y))


max.iter <- 250
out <- gbm.sc(Y1,oos.Y=Y2,M=20,max.iter=0)

U.init <- out$U %*% diag(sqrt(out$D), nrow=20)
V.init <- out$V %*% diag(sqrt(out$D), nrow=20)




## GLM-PCA Avagrad
library(glmpca)
max.iter <- 250

time.3 <- c(1:max.iter)
ll.3 <- c(1:max.iter)

fit <- glmpca(Y1,L=20,Y.oos=Y2,ctl=list(maxIter=max.iter),
              init=list(loadings=U.init, factors=V.init))
time.3 <- fit$mylist$time[-1]
ll.3 <- fit$mylist$LL

max.iter <- 100
fit <- glmpca(Y1,L=20,Y.oos=Y2,optimizer="fisher",ctl=list(verbose=TRUE,maxIter=max.iter),
              init=list(loadings=U.init, factors=V.init))
#No fisher for this one


time.4 <- fit$mylist$time[-1]
ll.4 <- fit$mylist$LL

set.seed(1)
max.iter <- 500
fit <- glmpca(Y1,Y.oos=Y2,L=20,minibatch="stochastic",ctl=list(verbose=TRUE,max.iter=max.iter,batch_size=400),
              init=list(loadings=U.init, factors=V.init))


time.5 <- fit$mylist$time[-1]
ll.5 <- fit$mylist$LL


Time <- list(time.3,time.4, time.5)
LL <- list(ll.3,ll.4, ll.5)


save(Time,file="../data/zhengmix_time_gbmInit.RData")
save(LL,file="../data/zhengmix_LL_gbmInit.RData")




















##BLISH

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
library(scGBM)
library(fastglm)
library(bigmemory)


#Sco <- readRDS("../../data/blish.RDS")
#Sco <- NormalizeData(Sco, assay="RNA")
#Sco <- FindVariableFeatures(Sco,assay="RNA",nfeatures=2000)
Y <- readRDS("../../data/blish_counts.RDS")
Y <- as.matrix(Y)

Y <- Y[rowSums(Y) >= 50,]
ds <- scGBM:::data.split(Y,p=0.5)
Y1 <- ds$Y1; Y2 <- ds$Y2
I <- nrow(Y); J <- ncol(Y)
print(dim(Y))

max.iter <- 250
out <- gbm.sc(Y1,oos.Y=Y2,M=20,max.iter=0)

U.init <- out$U %*% diag(sqrt(out$D), nrow=20)
V.init <- out$V %*% diag(sqrt(out$D), nrow=20)



## GLM-PCA Avagrad
library(glmpca)
max.iter <- 250

time.3 <- c(1:max.iter)
ll.3 <- c(1:max.iter)

fit <- glmpca(Y1,L=20,Y.oos=Y2,ctl=list(maxIter=max.iter),
              init=list(loadings=U.init, factors=V.init))
time.3 <- fit$mylist$time[-1]
ll.3 <- fit$mylist$LL

max.iter <- 100
fit <- glmpca(Y1,L=20,Y.oos=Y2,optimizer="fisher",ctl=list(verbose=TRUE,maxIter=max.iter),
              init=list(loadings=U.init, factors=V.init))
#No fisher for this one


time.4 <- fit$mylist$time[-1]
ll.4 <- fit$mylist$LL

set.seed(1)
max.iter <- 500
fit <- glmpca(Y1,Y.oos=Y2,L=20,minibatch="stochastic",ctl=list(verbose=TRUE,max.iter=max.iter,batch_size=400),
              init=list(loadings=U.init, factors=V.init))


time.5 <- fit$mylist$time[-1]
ll.5 <- fit$mylist$LL


Time <- list(time.3,time.4, time.5)
LL <- list(ll.3,ll.4, ll.5)


save(Time,file="../data/blish_time_gbmInit.RData")
save(LL,file="../data/blish_LL_gbmInit.RData")











