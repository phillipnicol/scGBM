

### Subsample
library(DuoClustering2018)
library(pdfCluster)
library(glmpca)
library(scGBM)
library(Seurat)
library(scGBM)

set.seed(42)

## 8eq
sce <- DuoClustering2018::sce_full_Zhengmix8eq()
phenoid <- sce$phenoid
Y <- sce@assays@data$counts
Y <- Y[rowSums(Y) >= 5,]


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



zhengmix_countsplit <- function(Y) {
  # Generate count split data
  #set.seed(1)
  ds <- scGBM:::data.split(Y, 0.5)

  Y1 <- ds$Y1
  Y2 <- ds$Y2

  Y1 <- Y1[rowSums(Y1) >= 5, ]
  Y2 <- Y2[rowSums(Y2) >= 5, ]

  # scGBM
  out <- gbm.sc(Y1, M = 20, max.iter = 250, tol = 10^-5, sigma = 10)
  rownames(Y1) <- 1:nrow(Y1)
  colnames(Y1) <- colnames(Y)
  Sco <- CreateSeuratObject(counts = Y1)
  colnames(out$scores) <- 1:20
  rownames(out$scores) <- colnames(Y1)
  Sco[["gbm"]] <- CreateDimReducObject(embeddings = out$scores, key = "GBM_")
  Sco <- FindNeighbors(Sco, reduction = "gbm")
  Sco <- FindClusters(Sco)
  cluster1 <- Sco$seurat_clusters[phenoid == "cd14.monocytes"]

  out <- gbm.sc(Y2, M = 20, max.iter = 250, tol = 10^-5, sigma = 10)
  rownames(Y2) <- 1:nrow(Y2)
  colnames(Y2) <- colnames(Y)
  Sco <- CreateSeuratObject(counts = Y2)
  colnames(out$scores) <- 1:20
  rownames(out$scores) <- colnames(Y2)
  Sco[["gbm"]] <- CreateDimReducObject(embeddings = out$scores, key = "GBM_")
  Sco <- FindNeighbors(Sco, reduction = "gbm")
  Sco <- FindClusters(Sco)
  cluster2 <- Sco$seurat_clusters[phenoid == "cd14.monocytes"]

  gbm.mono <- sum(choose(table(cluster1, cluster2), 2)) / choose(600, 2)

  # APR+PCA
  Sco <- CreateSeuratObject(counts = Y1)
  Sco <- SCTransform(Sco, vst.flavor = "v1", method = "offset", variable.features.n = nrow(Y1))
  Sco <- RunPCA(Sco, npcs = 20)
  Sco <- FindNeighbors(Sco)
  Sco <- FindClusters(Sco)
  cluster1 <- Sco$seurat_clusters[phenoid == "cd14.monocytes"]

  Sco <- CreateSeuratObject(counts = Y2)
  Sco <- SCTransform(Sco, vst.flavor = "v1", method = "offset", variable.features.n = nrow(Y2))
  Sco <- RunPCA(Sco, npcs = 20)
  Sco <- FindNeighbors(Sco)
  Sco <- FindClusters(Sco)
  cluster2 <- Sco$seurat_clusters[phenoid == "cd14.monocytes"]

  apr.mono <- sum(choose(table(cluster1, cluster2), 2)) / choose(600, 2)

  # LOG+Scale+PCA
  Sco <- CreateSeuratObject(counts = Y1)
  Sco <- NormalizeData(Sco)
  Sco <- FindVariableFeatures(Sco, nfeatures = nrow(Y1))
  Sco <- ScaleData(Sco)
  Sco <- RunPCA(Sco, npcs = 20)
  Sco <- FindNeighbors(Sco)
  Sco <- FindClusters(Sco)
  cluster1 <- Sco$seurat_clusters[phenoid == "cd14.monocytes"]

  Sco <- CreateSeuratObject(counts = Y2)
  Sco <- NormalizeData(Sco)
  Sco <- FindVariableFeatures(Sco, nfeatures = nrow(Y2))
  Sco <- ScaleData(Sco)
  Sco <- RunPCA(Sco, npcs = 20)
  Sco <- FindNeighbors(Sco)
  Sco <- FindClusters(Sco)
  cluster2 <- Sco$seurat_clusters[phenoid == "cd14.monocytes"]

  lpca.mono <- sum(choose(table(cluster1, cluster2), 2)) / choose(600, 2)

  # SCT+PCA
  Sco <- CreateSeuratObject(counts = Y1)
  Sco <- SCTransform(Sco, variable.features.n = nrow(Y1))
  Sco <- RunPCA(Sco, npcs = 20)
  Sco <- FindNeighbors(Sco)
  Sco <- FindClusters(Sco)
  cluster1 <- Sco$seurat_clusters[phenoid == "cd14.monocytes"]

  Sco <- CreateSeuratObject(counts = Y2)
  Sco <- SCTransform(Sco, variable.features.n = nrow(Y2))
  Sco <- RunPCA(Sco, npcs = 20)
  Sco <- FindNeighbors(Sco)
  Sco <- FindClusters(Sco)
  cluster2 <- Sco$seurat_clusters[phenoid == "cd14.monocytes"]

  sct.mono <- sum(choose(table(cluster1, cluster2), 2)) / choose(600, 2)

  # LOG+PCA
  Sco <- CreateSeuratObject(counts = Y1)
  L <- median(colSums(Y1))
  YL <- log(sweep(Y1, MARGIN = 2, STATS = L^-1 * colSums(Y1), FUN = "/") + 1)
  my.pca <- irlba::prcomp_irlba(t(YL), n = 20)
  rownames(my.pca$x) <- colnames(Y1)
  Sco[["lpca"]] <- CreateDimReducObject(embeddings = my.pca$x, key = "LPCA_")
  Sco <- FindNeighbors(Sco, reduction = "lpca")
  Sco <- FindClusters(Sco)
  cluster1 <- Sco$seurat_clusters[phenoid == "cd14.monocytes"]

  Sco <- CreateSeuratObject(counts = Y2)
  L <- median(colSums(Y2))
  YL <- log(sweep(Y2, MARGIN = 2, STATS = L^-1 * colSums(Y1), FUN = "/") + 1)
  my.pca <- irlba::prcomp_irlba(t(YL), n = 20)
  rownames(my.pca$x) <- colnames(Y2)
  Sco[["lpca"]] <- CreateDimReducObject(embeddings = my.pca$x, key = "LPCA_")
  Sco <- FindNeighbors(Sco, reduction = "lpca")
  Sco <- FindClusters(Sco)
  cluster2 <- Sco$seurat_clusters[phenoid == "cd14.monocytes"]

  lpca_ns.mono <- sum(choose(table(cluster1, cluster2), 2)) / choose(600, 2)

  # GLMPCA (SGD)
  set.seed(1)
  fit <- glmpca(Y1, L = 20, minibatch = "stochastic", ctl = list(batch_size = 400))
  Sco <- CreateSeuratObject(counts = Y1)
  Sco[["glmpca"]] <- CreateDimReducObject(embeddings = as.matrix(fit$res$factors), key = "GLMPCA_")
  Sco <- FindNeighbors(Sco, reduction = "glmpca")
  Sco <- FindClusters(Sco)
  cluster1 <- Sco$seurat_clusters[phenoid == "cd14.monocytes"]

  set.seed(1)
  fit <- glmpca(Y2, L = 20, minibatch = "stochastic", ctl = list(batch_size = 400))
  Sco <- CreateSeuratObject(counts = Y2)
  Sco[["glmpca"]] <- CreateDimReducObject(embeddings = as.matrix(fit$res$factors), key = "GLMPCA_")
  Sco <- FindNeighbors(Sco, reduction = "glmpca")
  Sco <- FindClusters(Sco)
  cluster2 <- Sco$seurat_clusters[phenoid == "cd14.monocytes"]

  glmpca.mono <- sum(choose(table(cluster1, cluster2), 2)) / choose(600, 2)

  # Return results
  return(c(
    gbm.mono = gbm.mono,
    apr.mono = apr.mono,
    lpca.mono = lpca.mono,
    sct.mono = sct.mono,
    lpca_ns.mono = lpca_ns.mono,
    glmpca.mono = glmpca.mono
  ))
}


