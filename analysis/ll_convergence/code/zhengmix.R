
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
out <- gbm.sc(Y1,oos.Y=Y2,M=20,max.iter=max.iter,tol=10^{-4},infer.beta=TRUE,time.by.iter = TRUE)
print(out$ll.oos)

gbm_umap <- umap::umap(out$scores)$layout
saveRDS(gbm_umap, file="../data/gbm_zhengmix_umap.RDS")

time.1 <- out$time
ll.1 <- out$ll.oos

saveRDS(out$scores, file="../data/gbm_zhengmix_embedding.RDS")

## GBM SC PAR
## GBM SC PAR
library(parallel)
library(fastglm)

niter <- length(seq(25,250,by=25))
proj.time.all <- matrix(0, nrow=niter,ncol=10) #ten reps
proj.ll.all <- proj.time.all

for(i in 1:10) {
  proj.time <- c()
  proj.ll <- c()
  #subset <- sample(1:ncol(Y),size=400,replace=FALSE)

  for(k in seq(25,250,by=25)) {
    print(k)
    start <- Sys.time()
    out <- gbm.sc(Y1,M=20,max.iter=k,subset=400,ncores=12,tol=10^{-4})
    end <- Sys.time()
    proj.time <- c(proj.time,difftime(end,start,units="secs")[[1]])

    Alpha <- matrix(out$alpha,nrow=I,ncol=J)
    Beta <- matrix(out$beta,nrow=I,ncol=J,byrow=TRUE)
    W <- exp(Alpha+Beta+out$U %*% t(out$scores))
    proj.ll <- c(proj.ll,sum(Y2*log(W) - W))
    print(sum(Y2*log(W) - W))
  }

  proj.time.all[,i] <- proj.time
  proj.ll.all[,i] <- proj.ll
}
ll.2 <- proj.ll.all
time.2 <- proj.time.all
print(ll.2)

saveRDS(out$scores, file="../data/gbm_proj_zhengmix_embedding.RDS")




## GLM-PCA Avagrad
library(glmpca)
max.iter <- 250

time.3 <- c(1:max.iter)
ll.3 <- c(1:max.iter)

fit <- glmpca(Y1,L=20,Y.oos=Y2,ctl=list(maxIter=max.iter))
time.3 <- fit$mylist$time[-1]
ll.3 <- fit$mylist$LL

saveRDS(fit$res$factors, file="../data/glmpca_avagrad_zhengmix_embedding.RDS")

max.iter <- 100
fit <- glmpca(Y1,L=20,Y.oos=Y2,optimizer="fisher",ctl=list(verbose=TRUE,maxIter=max.iter))

saveRDS(fit$res$factors, file="../data/glmpca_fisher_zhengmix_embedding.RDS")

time.4 <- fit$mylist$time[-1]
ll.4 <- fit$mylist$LL

set.seed(1)
max.iter <- 500
fit <- glmpca(Y1,Y.oos=Y2,L=20,minibatch="stochastic",ctl=list(verbose=TRUE,max.iter=max.iter,batch_size=400))

saveRDS(fit$res$factors, file="../data/glmpca_sgd_zhengmix_embedding.RDS")

time.5 <- fit$mylist$time[-1]
ll.5 <- fit$mylist$LL


nalgos <- 3
Time <- list(time.1,time.2,time.3,time.4, time.5)
LL <- list(ll.1,ll.2,ll.3,ll.4, ll.5)


save(Time,file="../data/zhengmix_time.RData")
save(LL,file="../data/zhengmix_LL.RData")


### Comparison to standard pipeline

Y <- Y1 #Use the one that was compared against

## Comparison to LOG+PCA

L <- median(colSums(Y))
YL <- log(sweep(Y,MARGIN=2,STATS=L^{-1}*colSums(Y),FUN="/") + 1)
my.pca <- irlba::prcomp_irlba(t(YL),n=20)
saveRDS(my.pca$x, "../data/logpca_zhemgmix_embedding.RDS")

pca.umap <- umap::umap(my.pca$x)$layout
saveRDS(pca.umap, "../data/logpca_zhengmix_umap.RDS")

## Comparison to APR

apr <- sctransform::vst(Y, method="offset")
pca.apr <- irlba::prcomp_irlba(t(apr$y),n=20)
saveRDS(pca.apr$x, "../data/apr_zhengmix_embedding.RDS")

apr.umap <- umap::umap(pca.apr$x)$layout
saveRDS(apr.umap, "../data/apr_zhengmix_umap.RDS")

## Comparison to SCT
Sco <- CreateSeuratObject(counts=Y)
Sco <- SCTransform(Sco)
Sco <- RunPCA(Sco)
sct <- Sco@reductions$pca@cell.embeddings
saveRDS(sct, "../data/sct_zhengmix_embedding.RDS")

sct.umap <- umap::umap(sct)$layout
saveRDS(sct.umap, "../data/sct_zhengmix_umap.RDS")

## Comparison to Log+Scale+PCA

Sco <- CreateSeuratObject(counts=Y)
Sco <- NormalizeData(Sco)
Sco <- FindVariableFeatures(Sco)
Sco <- ScaleData(Sco)
Sco <- RunPCA(Sco)

lpca <- Sco@reductions$pca@cell.embeddings
saveRDS(lpca, "../data/seurat_zhengmix_embedding.RDS")
lpca.scale.umap <- umap::umap(lpca)$layout
saveRDS(lpca.scale.umap, "../data/seurat_zhengmix_umap.RDS")



