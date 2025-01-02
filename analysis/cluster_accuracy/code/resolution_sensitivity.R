

run_clustering_resolution <- function(resolution) {

  ### Subsample
  library(DuoClustering2018)
  library(pdfCluster)
  library(glmpca)
  library(scGBM)
  library(Seurat)
  library(scGBM)
  library(ggplot2)

  set.seed(42)

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
  Sco <- FindClusters(Sco, resolution = resolution)
  gbm <- adj.rand.index(phenoid, Sco$seurat_clusters)



  fit <- glmpca(Y, L=20)
  Sco <- CreateSeuratObject(counts=Y)
  Sco[["glmpca"]] <- CreateDimReducObject(embeddings=as.matrix(fit$res$factors),key="GLMPCA_")
  Sco <- FindNeighbors(Sco,reduction = "glmpca")
  Sco <- FindClusters(Sco, resolution = resolution)
  glmpca <- adj.rand.index(phenoid, Sco$seurat_clusters)


  Sco <- CreateSeuratObject(counts=Y)
  Sco <- NormalizeData(Sco)
  Sco <- FindVariableFeatures(Sco,nfeatures=nrow(Y))
  Sco <- ScaleData(Sco)
  Sco <- RunPCA(Sco,npcs=20)
  Sco <- FindNeighbors(Sco)
  Sco <- FindClusters(Sco, resolution = resolution)
  l2pca <- adj.rand.index(phenoid, Sco$seurat_clusters)


  Sco <- CreateSeuratObject(counts=Y)
  Sco <- SCTransform(Sco, variable.features.n = nrow(Y))
  Sco <- RunPCA(Sco,npcs=20)
  Sco <- FindNeighbors(Sco)
  Sco <- FindClusters(Sco, resolution = resolution)
  sct <- adj.rand.index(phenoid, Sco$seurat_clusters)

  Sco <- CreateSeuratObject(counts=Y)
  Sco <- SCTransform(Sco, vst.flavor = "v1", method="offset",variable.features.n = nrow(Y))
  Sco <- RunPCA(Sco,npcs=20)
  Sco <- FindNeighbors(Sco)
  Sco <- FindClusters(Sco, resolution = resolution)
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
    Sco <- FindClusters(Sco, resolution = resolution)
    proj_res[j] <- adj.rand.index(phenoid, Sco$seurat_clusters)
  }
  gbmproj <- mean(proj_res)

  ### LOG-PCA
  Sco <- CreateSeuratObject(counts=Y)
  L <- median(colSums(Y))
  YL <- log(sweep(Y,MARGIN=2,STATS=L^{-1}*colSums(Y),FUN="/") + 1)
  my.pca <- irlba::prcomp_irlba(t(YL), n=20)
  rownames(my.pca$x) <- colnames(Y)
  Sco[["lpca"]] <- CreateDimReducObject(embeddings=my.pca$x,key="LPCA_")
  Sco <- FindNeighbors(Sco,reduction = "lpca")
  Sco <- FindClusters(Sco, resolution = resolution)
  lpca_noscale <- adj.rand.index(phenoid, Sco$seurat_clusters)

  results <- c(gbm,gbmproj, l2pca, sct,apr,glmpca,lpca_noscale)
  names(results) <- c("scGBM-full", "scGBM-proj", "Log+Scale+PCA", "SCT", "APR","GLM-PCA (SGD)",
                      "Log+PCA")

  saveRDS(results, file="../data/zhengmix8eq.RDS")



  ##### CLASS SIZE IMBALANCE #####

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
  Sco <- FindClusters(Sco, resolution = resolution)
  gbm <- adj.rand.index(phenoid, Sco$seurat_clusters)


  set.seed(1)
  fit <- glmpca(Y, L=20, minibatch="stochastic", ctl=list(batch_size=400))
  Sco <- CreateSeuratObject(counts=Y)
  Sco[["glmpca"]] <- CreateDimReducObject(embeddings=as.matrix(fit$res$factors),key="GLMPCA_")
  Sco <- FindNeighbors(Sco,reduction = "glmpca")
  Sco <- FindClusters(Sco, resolution = resolution)
  glmpca <- adj.rand.index(phenoid, Sco$seurat_clusters)

  ##Log+Scale+PCA
  Sco <- CreateSeuratObject(counts=Y)
  Sco <- NormalizeData(Sco)
  Sco <- FindVariableFeatures(Sco,nfeatures=nrow(Y))
  Sco <- ScaleData(Sco)
  Sco <- RunPCA(Sco,npcs=20)
  Sco <- FindNeighbors(Sco)
  Sco <- FindClusters(Sco, resolution = resolution)
  l2pca <- adj.rand.index(phenoid, Sco$seurat_clusters)

  ##SCTRANSFORM
  Sco <- CreateSeuratObject(counts=Y)
  Sco <- SCTransform(Sco, variable.features.n = nrow(Y))
  Sco <- RunPCA(Sco,npcs=20)
  Sco <- FindNeighbors(Sco)
  Sco <- FindClusters(Sco, resolution = resolution)
  sct <- adj.rand.index(phenoid, Sco$seurat_clusters)


  ###APR
  Sco <- CreateSeuratObject(counts=Y)
  Sco <- SCTransform(Sco, vst.flavor = "v1", method="offset",variable.features.n = nrow(Y))
  Sco <- RunPCA(Sco,npcs=20)
  Sco <- FindNeighbors(Sco)
  Sco <- FindClusters(Sco, resolution = resolution)
  #apr <- adj.rand.index(sce$phenoid, Sco$seurat_clusters)
  apr <- adj.rand.index(phenoid, Sco$seurat_clusters)

  ##SCTRANSFORM
  Sco <- CreateSeuratObject(counts=Y)
  Sco <- SCTransform(Sco, variable.features.n = nrow(Y))
  Sco <- RunPCA(Sco,npcs=20)
  Sco <- FindNeighbors(Sco)
  Sco <- FindClusters(Sco, resolution = resolution)
  sct <- adj.rand.index(phenoid, Sco$seurat_clusters)


  ### LOG-PCA
  Sco <- CreateSeuratObject(counts=Y)
  L <- median(colSums(Y))
  YL <- log(sweep(Y,MARGIN=2,STATS=L^{-1}*colSums(Y),FUN="/") + 1)
  my.pca <- irlba::prcomp_irlba(t(YL), n=20)
  rownames(my.pca$x) <- colnames(Y)
  Sco[["lpca"]] <- CreateDimReducObject(embeddings=my.pca$x,key="LPCA_")
  Sco <- FindNeighbors(Sco,reduction = "lpca")
  Sco <- FindClusters(Sco, resolution = resolution)
  lpca_noscale <- adj.rand.index(phenoid, Sco$seurat_clusters)

  library(fastglm)
  proj_res <- rep(0, 10)
  for(j in 1:10) {
    outproj <- gbm.sc(Y,M=20,subset=400,ncores=8)
    Sco <- CreateSeuratObject(counts=Y)
    colnames(outproj$scores) <- 1:20
    rownames(outproj$scores) <- colnames(Y)
    Sco[["gbm"]] <- CreateDimReducObject(embeddings=outproj$scores,key="GBM_")
    Sco <- FindNeighbors(Sco,reduction = "gbm")
    Sco <- FindClusters(Sco, resolution = resolution)
    proj_res[j] <- adj.rand.index(phenoid, Sco$seurat_clusters)
  }
  gbmproj <- mean(proj_res)

  results_subsampled <- c(gbm,gbmproj, l2pca, sct,apr,glmpca,lpca_noscale)
  names(results_subsampled) <- c("scGBM-full", "scGBM-proj", "Log+Scale+PCA", "SCT", "APR","GLM-PCA (SGD)",
                                 "Log+PCA")

  res <- list()
  res$results <- results
  res$results_subsampled <- results_subsampled
  return(res)
}


resolutions <- seq(0.3, 1.5, by=0.1)

res_mat <- matrix(0, nrow=length(resolutions), ncol=7)
colnames(res_mat) <- c("scGBM-full", "scGBM-proj", "Log+Scale+PCA", "SCT", "APR","GLM-PCA (SGD)",
                       "Log+PCA")

res_subsampled_mat <- matrix(0, nrow=length(resolutions), ncol=7)
colnames(res_subsampled_mat) <- c("scGBM-full", "scGBM-proj", "Log+Scale+PCA", "SCT", "APR","GLM-PCA (SGD)",
                       "Log+PCA")


for(i in 1:length(resolutions)) {
  res_list <- run_clustering_resolution(resolution=resolutions[i])

  res_mat[i,] <- res_list$results
  res_subsampled_mat[i,] <- res_list$results_subsampled

  cat("RESULTS: ", res_mat[i,], " ", res_subsampled_mat[i,], "\n")
}

saveRDS(res_mat, file="../data/resolution_sensitivity_full.RDS")
saveRDS(res_subsampled_mat, file="../data/resolution_sensitivity_subsampled.RDS")
