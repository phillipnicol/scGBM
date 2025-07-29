
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
out <- gbm.sc(Y1,oos.Y=Y2,M=20,max.iter=max.iter,tol=0,
              infer.beta=TRUE,
              time.by.iter = TRUE,
              factor.init = "near-zero")
print(out$ll.oos)

time.1 <- out$time
ll.1 <- out$ll.oos

##ZHENGMIX


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

max.iter <- 100
out <- gbm.sc(Y1,oos.Y=Y2,M=20,max.iter=max.iter,tol=10^{-4},infer.beta=TRUE,
              time.by.iter = TRUE,
              factor.init = "near-zero")
print(out$ll.oos)


time.2 <- out$time
ll.2 <- out$ll.oos


near.zero.results <- list(time.zhengmix = time.2,
                          ll.zhengmix = ll.2,
                          time.blish = time.1,
                          ll.blish = ll.2)

save(near.zero.results,
     file = "../data/near_zero_initialization.RData")


