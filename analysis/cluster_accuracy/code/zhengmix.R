
### Subsample
library(DuoClustering2018)
library(pdfCluster)
library(glmpca)
library(scGBM)
library(Seurat)
library(scGBM)

## 8eq
sce <- DuoClustering2018::sce_full_Zhengmix8eq()
phenoid <- sce$phenoid
Y <- sce@assays@data$counts
Y <- Y[rowSums(Y) >= 5,]

out <- gbm.sc(Y, M=20, max.iter=250, tol=10^{-5},sigma=10)

Sco <- CreateSeuratObject(counts=Y)
colnames(out$scores) <- 1:20
Sco[["gbm"]] <- CreateDimReducObject(embeddings=out$scores,key="GBM_")
Sco <- FindNeighbors(Sco,reduction = "gbm")
Sco <- FindClusters(Sco)
gbm <- adj.rand.index(phenoid, Sco$seurat_clusters)

fit <- glmpca(Y, L=20)
Sco <- CreateSeuratObject(counts=Y)
Sco[["glmpca"]] <- CreateDimReducObject(embeddings=as.matrix(fit$res$factors),key="GLMPCA_")
Sco <- FindNeighbors(Sco,reduction = "glmpca")
Sco <- FindClusters(Sco)
glmpca <- adj.rand.index(phenoid, Sco$seurat_clusters)


Sco <- CreateSeuratObject(counts=Y)
Sco <- NormalizeData(Sco)
Sco <- FindVariableFeatures(Sco,nfeatures=nrow(Y))
Sco <- ScaleData(Sco)
Sco <- RunPCA(Sco,npcs=20)
Sco <- FindNeighbors(Sco)
Sco <- FindClusters(Sco)
l2pca <- adj.rand.index(phenoid, Sco$seurat_clusters)


Sco <- CreateSeuratObject(counts=Y)
Sco <- SCTransform(Sco, variable.features.n = nrow(Y))
Sco <- RunPCA(Sco,npcs=20)
Sco <- FindNeighbors(Sco)
Sco <- FindClusters(Sco)
sct <- adj.rand.index(phenoid, Sco$seurat_clusters)

Sco <- CreateSeuratObject(counts=Y)
Sco <- SCTransform(Sco, vst.flavor = "v1", method="offset",variable.features.n = nrow(Y))
Sco <- RunPCA(Sco,npcs=20)
Sco <- FindNeighbors(Sco)
Sco <- FindClusters(Sco)
#apr <- adj.rand.index(sce$phenoid, Sco$seurat_clusters)
apr <- adj.rand.index(phenoid, Sco$seurat_clusters)

library(fastglm)
proj_res <- rep(0, 10)
for(j in 1:10) {
  outproj <- gbm.sc(Y,M=20,subset=640,ncores=8)
  Sco <- CreateSeuratObject(counts=Y)
  colnames(outproj$scores) <- 1:20
  rownames(outproj$scores) <- colnames(Y)
  Sco[["gbm"]] <- CreateDimReducObject(embeddings=outproj$scores,key="GBM_")
  Sco <- FindNeighbors(Sco,reduction = "gbm")
  Sco <- FindClusters(Sco)
  proj_res[j] <- adj.rand.index(phenoid, Sco$seurat_clusters)
}
gbmproj <- mean(proj_res)

results <- c(gbm,gbmproj, l2pca, sct,apr,glmpca)
names(results) <- c("scGBM-full", "scGBM-proj", "log+scale+PCA", "SCT", "APR","GLM-PCA")

saveRDS(results, file="../data/zhengmix8eq.RDS")



### 8uneq
sizes <- c(100, 100, 600, 100, 100, 100, 100, 100)
set.seed(1)
ixs <- c()
phenoid <- c()
for(i in 1:length(unique(sce$phenoid))) {
  ixs <- c(ixs, sample(which(sce$phenoid == unique(sce$phenoid)[i]), size=sizes[i]))
  phenoid <- c(phenoid, rep(unique(sce$phenoid)[i], sizes[i]))
}

Y <- Y[,ixs]
Y <- Y[rowSums(Y) >= 5,]



out <- gbm.sc(Y, M=20, max.iter=250, tol=10^{-5},sigma=10)

Sco <- CreateSeuratObject(counts=Y)
colnames(out$scores) <- 1:20
Sco[["gbm"]] <- CreateDimReducObject(embeddings=out$scores,key="GBM_")
Sco <- FindNeighbors(Sco,reduction = "gbm")
Sco <- FindClusters(Sco)
gbm <- adj.rand.index(phenoid, Sco$seurat_clusters)

fit <- glmpca(Y, L=20)
Sco <- CreateSeuratObject(counts=Y)
Sco[["glmpca"]] <- CreateDimReducObject(embeddings=as.matrix(fit$res$factors),key="GLMPCA_")
Sco <- FindNeighbors(Sco,reduction = "glmpca")
Sco <- FindClusters(Sco)
glmpca <- adj.rand.index(phenoid, Sco$seurat_clusters)


Sco <- CreateSeuratObject(counts=Y)
Sco <- NormalizeData(Sco)
Sco <- FindVariableFeatures(Sco,nfeatures=nrow(Y))
Sco <- ScaleData(Sco)
Sco <- RunPCA(Sco,npcs=20)
Sco <- FindNeighbors(Sco)
Sco <- FindClusters(Sco)
l2pca <- adj.rand.index(phenoid, Sco$seurat_clusters)


Sco <- CreateSeuratObject(counts=Y)
Sco <- SCTransform(Sco, variable.features.n = nrow(Y))
Sco <- RunPCA(Sco,npcs=20)
Sco <- FindNeighbors(Sco)
Sco <- FindClusters(Sco)
sct <- adj.rand.index(phenoid, Sco$seurat_clusters)

Sco <- CreateSeuratObject(counts=Y)
Sco <- SCTransform(Sco, vst.flavor = "v1", method="offset",variable.features.n = nrow(Y))
Sco <- RunPCA(Sco,npcs=20)
Sco <- FindNeighbors(Sco)
Sco <- FindClusters(Sco)
#apr <- adj.rand.index(sce$phenoid, Sco$seurat_clusters)
apr <- adj.rand.index(phenoid, Sco$seurat_clusters)

library(fastglm)
proj_res <- rep(0, 10)
for(j in 1:10) {
  outproj <- gbm.sc(Y,M=20,subset=400,ncores=8)
  Sco <- CreateSeuratObject(counts=Y)
  colnames(outproj$scores) <- 1:20
  rownames(outproj$scores) <- colnames(Y)
  Sco[["gbm"]] <- CreateDimReducObject(embeddings=outproj$scores,key="GBM_")
  Sco <- FindNeighbors(Sco,reduction = "gbm")
  Sco <- FindClusters(Sco)
  proj_res[j] <- adj.rand.index(phenoid, Sco$seurat_clusters)
}
gbmproj <- mean(proj_res)

results <- c(gbm,gbmproj, l2pca, sct,apr,glmpca)
names(results) <- c("scGBM-full", "scGBM-proj", "log+scale+PCA", "SCT", "APR","GLM-PCA")

saveRDS(results, file="../data/zhengmix8uneq.RDS")





## 4eq
sce <- DuoClustering2018::sce_full_Zhengmix4eq()
phenoid <- sce$phenoid
Y <- sce@assays@data$counts
Y <- Y[rowSums(Y) >= 5,]

out <- gbm.sc(Y, M=20, max.iter=250, tol=10^{-5},sigma=10)

Sco <- CreateSeuratObject(counts=Y)
colnames(out$scores) <- 1:20
Sco[["gbm"]] <- CreateDimReducObject(embeddings=out$scores,key="GBM_")
Sco <- FindNeighbors(Sco,reduction = "gbm")
Sco <- FindClusters(Sco)
gbm <- adj.rand.index(phenoid, Sco$seurat_clusters)

fit <- glmpca(Y, L=20)
Sco <- CreateSeuratObject(counts=Y)
Sco[["glmpca"]] <- CreateDimReducObject(embeddings=as.matrix(fit$res$factors),key="GLMPCA_")
Sco <- FindNeighbors(Sco,reduction = "glmpca")
Sco <- FindClusters(Sco)
glmpca <- adj.rand.index(phenoid, Sco$seurat_clusters)


Sco <- CreateSeuratObject(counts=Y)
Sco <- NormalizeData(Sco)
Sco <- FindVariableFeatures(Sco,nfeatures=nrow(Y))
Sco <- ScaleData(Sco)
Sco <- RunPCA(Sco,npcs=20)
Sco <- FindNeighbors(Sco)
Sco <- FindClusters(Sco)
l2pca <- adj.rand.index(phenoid, Sco$seurat_clusters)


Sco <- CreateSeuratObject(counts=Y)
Sco <- SCTransform(Sco, variable.features.n = nrow(Y))
Sco <- RunPCA(Sco,npcs=20)
Sco <- FindNeighbors(Sco)
Sco <- FindClusters(Sco)
sct <- adj.rand.index(phenoid, Sco$seurat_clusters)

Sco <- CreateSeuratObject(counts=Y)
Sco <- SCTransform(Sco, vst.flavor = "v1", method="offset",variable.features.n = nrow(Y))
Sco <- RunPCA(Sco,npcs=20)
Sco <- FindNeighbors(Sco)
Sco <- FindClusters(Sco)
#apr <- adj.rand.index(sce$phenoid, Sco$seurat_clusters)
apr <- adj.rand.index(phenoid, Sco$seurat_clusters)

library(fastglm)
proj_res <- rep(0, 10)
for(j in 1:10) {
  outproj <- gbm.sc(Y,M=20,subset=400,ncores=8)
  Sco <- CreateSeuratObject(counts=Y)
  colnames(outproj$scores) <- 1:20
  rownames(outproj$scores) <- colnames(Y)
  Sco[["gbm"]] <- CreateDimReducObject(embeddings=outproj$scores,key="GBM_")
  Sco <- FindNeighbors(Sco,reduction = "gbm")
  Sco <- FindClusters(Sco)
  proj_res[j] <- adj.rand.index(phenoid, Sco$seurat_clusters)
}
gbmproj <- mean(proj_res)

results <- c(gbm,gbmproj, l2pca, sct,apr,glmpca)
names(results) <- c("scGBM-full", "scGBM-proj", "log+scale+PCA", "SCT", "APR","GLM-PCA")

saveRDS(results, file="../data/zhengmix4eq.RDS")


## 4uneq

sce <- DuoClustering2018::sce_filteredExpr10_Zhengmix4uneq()
phenoid <- sce$phenoid
Y <- sce@assays@data$counts
Y <- Y[rowSums(Y) >= 5,]

out <- gbm.sc(Y, M=20, max.iter=250, tol=10^{-5},sigma=10)

Sco <- CreateSeuratObject(counts=Y)
colnames(out$scores) <- 1:20
Sco[["gbm"]] <- CreateDimReducObject(embeddings=out$scores,key="GBM_")
Sco <- FindNeighbors(Sco,reduction = "gbm")
Sco <- FindClusters(Sco)
gbm <- adj.rand.index(phenoid, Sco$seurat_clusters)

fit <- glmpca(Y, L=20)
Sco <- CreateSeuratObject(counts=Y)
Sco[["glmpca"]] <- CreateDimReducObject(embeddings=as.matrix(fit$res$factors),key="GLMPCA_")
Sco <- FindNeighbors(Sco,reduction = "glmpca")
Sco <- FindClusters(Sco)
glmpca <- adj.rand.index(phenoid, Sco$seurat_clusters)


Sco <- CreateSeuratObject(counts=Y)
Sco <- NormalizeData(Sco)
Sco <- FindVariableFeatures(Sco,nfeatures=nrow(Y))
Sco <- ScaleData(Sco)
Sco <- RunPCA(Sco,npcs=20)
Sco <- FindNeighbors(Sco)
Sco <- FindClusters(Sco)
l2pca <- adj.rand.index(phenoid, Sco$seurat_clusters)


Sco <- CreateSeuratObject(counts=Y)
Sco <- SCTransform(Sco, variable.features.n = nrow(Y))
Sco <- RunPCA(Sco,npcs=20)
Sco <- FindNeighbors(Sco)
Sco <- FindClusters(Sco)
sct <- adj.rand.index(phenoid, Sco$seurat_clusters)

Sco <- CreateSeuratObject(counts=Y)
Sco <- SCTransform(Sco, vst.flavor = "v1", method="offset",variable.features.n = nrow(Y))
Sco <- RunPCA(Sco,npcs=20)
Sco <- FindNeighbors(Sco)
Sco <- FindClusters(Sco)
#apr <- adj.rand.index(sce$phenoid, Sco$seurat_clusters)
apr <- adj.rand.index(phenoid, Sco$seurat_clusters)

library(fastglm)
proj_res <- rep(0, 10)
for(j in 1:10) {
  outproj <- gbm.sc(Y,M=20,subset=640,ncores=8)
  Sco <- CreateSeuratObject(counts=Y)
  colnames(outproj$scores) <- 1:20
  rownames(outproj$scores) <- colnames(Y)
  Sco[["gbm"]] <- CreateDimReducObject(embeddings=outproj$scores,key="GBM_")
  Sco <- FindNeighbors(Sco,reduction = "gbm")
  Sco <- FindClusters(Sco)
  proj_res[j] <- adj.rand.index(phenoid, Sco$seurat_clusters)
}
gbmproj <- mean(proj_res)

results <- c(gbm,gbmproj, l2pca, sct,apr,glmpca)
names(results) <- c("scGBM-full", "scGBM-proj", "log+scale+PCA", "SCT", "APR","GLM-PCA")

saveRDS(results, file="../data/zhengmix4uneq.RDS")


