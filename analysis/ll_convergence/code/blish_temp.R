library(Seurat)
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


rownames(Y1) <- rownames(Y); colnames(Y1) <- colnames(Y)
Y <- Y1 #Use the one that was compared against


L <- median(colSums(Y))
YL <- log(sweep(Y,MARGIN=2,STATS=L^{-1}*colSums(Y),FUN="/") + 1)
my.pca <- irlba::prcomp_irlba(t(YL),n=20)
saveRDS(my.pca$x, "../data/logpca_blish_embedding.RDS")

pca.umap <- umap::umap(my.pca$x)$layout
saveRDS(pca.umap, "../data/logpca_blish_umap.RDS")

## Comparison to APR

apr <- sctransform::vst(Y, method="offset")
pca.apr <- irlba::prcomp_irlba(t(apr$y),n=20)
saveRDS(pca.apr$x, "../data/apr_blish_embedding.RDS")

apr.umap <- umap::umap(pca.apr$x)$layout
saveRDS(apr.umap, "../data/apr_blish_umap.RDS")

## Comparison to SCT
Sco <- CreateSeuratObject(counts=Y)
Sco <- SCTransform(Sco)
Sco <- RunPCA(Sco)
sct <- Sco@reductions$pca@cell.embeddings
saveRDS(sct, "../data/sct_blish_embedding.RDS")

sct.umap <- umap::umap(sct)$layout
saveRDS(sct.umap, "../data/sct_blish_umap.RDS")

## Comparison to Log+Scale+PCA

Sco <- CreateSeuratObject(counts=Y)
Sco <- NormalizeData(Sco)
Sco <- FindVariableFeatures(Sco)
Sco <- ScaleData(Sco)
Sco <- RunPCA(Sco)

lpca <- Sco@reductions$pca@cell.embeddings
saveRDS(lpca, "../data/seurat_blish_embedding.RDS")
lpca.scale.umap <- umap::umap(lpca)$layout
saveRDS(lpca.scale.umap, "../data/seurat_blish_umap.RDS")




