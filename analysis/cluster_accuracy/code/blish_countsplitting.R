

### Subsample
library(pdfCluster)
library(glmpca)
library(scGBM)
library(Seurat)
library(scGBM)

set.seed(42)



#Sco <- readRDS("../../data/blish.RDS")
#Sco <- NormalizeData(Sco, assay="RNA")
#Sco <- FindVariableFeatures(Sco,assay="RNA",nfeatures=2000)
Y <- readRDS("../../data/blish_counts.RDS")
Y <- as.matrix(Y)


blish_countsplit <- function(Y) {
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
  cluster1 <- Sco$seurat_clusters

  out <- gbm.sc(Y2, M = 20, max.iter = 250, tol = 10^-5, sigma = 10)
  rownames(Y2) <- 1:nrow(Y2)
  colnames(Y2) <- colnames(Y)
  Sco <- CreateSeuratObject(counts = Y2)
  colnames(out$scores) <- 1:20
  rownames(out$scores) <- colnames(Y2)
  Sco[["gbm"]] <- CreateDimReducObject(embeddings = out$scores, key = "GBM_")
  Sco <- FindNeighbors(Sco, reduction = "gbm")
  Sco <- FindClusters(Sco)
  cluster2 <- Sco$seurat_clusters

  gbm.mono <- sum(choose(table(cluster1, cluster2), 2)) / choose(ncol(Y),2)

  # APR+PCA
  Sco <- CreateSeuratObject(counts = Y1)
  Sco <- SCTransform(Sco, vst.flavor = "v1", method = "offset", variable.features.n = nrow(Y1))
  Sco <- RunPCA(Sco, npcs = 20)
  Sco <- FindNeighbors(Sco)
  Sco <- FindClusters(Sco)
  cluster1 <- Sco$seurat_clusters

  Sco <- CreateSeuratObject(counts = Y2)
  Sco <- SCTransform(Sco, vst.flavor = "v1", method = "offset", variable.features.n = nrow(Y2))
  Sco <- RunPCA(Sco, npcs = 20)
  Sco <- FindNeighbors(Sco)
  Sco <- FindClusters(Sco)
  cluster2 <- Sco$seurat_clusters

  apr.mono <- sum(choose(table(cluster1, cluster2), 2)) / choose(ncol(Y),2)

  # LOG+Scale+PCA
  Sco <- CreateSeuratObject(counts = Y1)
  Sco <- NormalizeData(Sco)
  Sco <- FindVariableFeatures(Sco, nfeatures = nrow(Y1))
  Sco <- ScaleData(Sco)
  Sco <- RunPCA(Sco, npcs = 20)
  Sco <- FindNeighbors(Sco)
  Sco <- FindClusters(Sco)
  cluster1 <- Sco$seurat_clusters

  Sco <- CreateSeuratObject(counts = Y2)
  Sco <- NormalizeData(Sco)
  Sco <- FindVariableFeatures(Sco, nfeatures = nrow(Y2))
  Sco <- ScaleData(Sco)
  Sco <- RunPCA(Sco, npcs = 20)
  Sco <- FindNeighbors(Sco)
  Sco <- FindClusters(Sco)
  cluster2 <- Sco$seurat_clusters

  lpca.mono <- sum(choose(table(cluster1, cluster2), 2)) / choose(ncol(Y),2)

  # SCT+PCA
  Sco <- CreateSeuratObject(counts = Y1)
  Sco <- SCTransform(Sco, variable.features.n = nrow(Y1))
  Sco <- RunPCA(Sco, npcs = 20)
  Sco <- FindNeighbors(Sco)
  Sco <- FindClusters(Sco)
  cluster1 <- Sco$seurat_clusters

  Sco <- CreateSeuratObject(counts = Y2)
  Sco <- SCTransform(Sco, variable.features.n = nrow(Y2))
  Sco <- RunPCA(Sco, npcs = 20)
  Sco <- FindNeighbors(Sco)
  Sco <- FindClusters(Sco)
  cluster2 <- Sco$seurat_clusters

  sct.mono <- sum(choose(table(cluster1, cluster2), 2)) / choose(ncol(Y),2)

  # LOG+PCA
  Sco <- CreateSeuratObject(counts = Y1)
  L <- median(colSums(Y1))
  YL <- log(sweep(Y1, MARGIN = 2, STATS = L^-1 * colSums(Y1), FUN = "/") + 1)
  my.pca <- irlba::prcomp_irlba(t(YL), n = 20)
  rownames(my.pca$x) <- colnames(Y1)
  Sco[["lpca"]] <- CreateDimReducObject(embeddings = my.pca$x, key = "LPCA_")
  Sco <- FindNeighbors(Sco, reduction = "lpca")
  Sco <- FindClusters(Sco)
  cluster1 <- Sco$seurat_clusters

  Sco <- CreateSeuratObject(counts = Y2)
  L <- median(colSums(Y2))
  YL <- log(sweep(Y2, MARGIN = 2, STATS = L^-1 * colSums(Y1), FUN = "/") + 1)
  my.pca <- irlba::prcomp_irlba(t(YL), n = 20)
  rownames(my.pca$x) <- colnames(Y2)
  Sco[["lpca"]] <- CreateDimReducObject(embeddings = my.pca$x, key = "LPCA_")
  Sco <- FindNeighbors(Sco, reduction = "lpca")
  Sco <- FindClusters(Sco)
  cluster2 <- Sco$seurat_clusters

  lpca_ns.mono <- sum(choose(table(cluster1, cluster2), 2)) / choose(ncol(Y),2)

  # GLMPCA (SGD)
  set.seed(1)
  fit <- glmpca(Y1, L = 20, minibatch = "stochastic", ctl = list(batch_size = 400))
  Sco <- CreateSeuratObject(counts = Y1)
  Sco[["glmpca"]] <- CreateDimReducObject(embeddings = as.matrix(fit$res$factors), key = "GLMPCA_")
  Sco <- FindNeighbors(Sco, reduction = "glmpca")
  Sco <- FindClusters(Sco)
  cluster1 <- Sco$seurat_clusters

  set.seed(1)
  fit <- glmpca(Y2, L = 20, minibatch = "stochastic", ctl = list(batch_size = 400))
  Sco <- CreateSeuratObject(counts = Y2)
  Sco[["glmpca"]] <- CreateDimReducObject(embeddings = as.matrix(fit$res$factors), key = "GLMPCA_")
  Sco <- FindNeighbors(Sco, reduction = "glmpca")
  Sco <- FindClusters(Sco)
  cluster2 <- Sco$seurat_clusters

  glmpca.mono <- sum(choose(table(cluster1, cluster2), 2)) / choose(ncol(Y),2)

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

iters <- 5
res <- matrix(0, nrow=5,ncol=6); colnames(res) <- c("gbm.mono",
                                                    "apr.mono",
                                                    "lpca.mono",
                                                    "sct.mono",
                                                    "lpca_ns.mono",
                                                    "glmpca.mono")


for(i in 1:5) {
  res[i, ] <- blish_countsplit(Y)
  cat("RESULTS: ", res[i,], "\n")
}

saveRDS(res, file="../data/blish_countsplit.RDS")

