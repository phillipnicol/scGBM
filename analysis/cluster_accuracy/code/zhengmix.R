setwd(here::here("analysis/cluster_accuracy/code"))

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

out <- gbm.sc(Y, M=20, max.iter=250, tol=10^{-5},sigma=10)

Sco <- CreateSeuratObject(counts=Y)
colnames(out$scores) <- 1:20
Sco[["gbm"]] <- CreateDimReducObject(embeddings=out$scores,key="GBM_")
Sco <- FindNeighbors(Sco,reduction = "gbm")
Sco <- FindClusters(Sco)
gbm <- adj.rand.index(phenoid, Sco$seurat_clusters)



set.seed(1)
fit <- glmpca(Y |> as.matrix(), L=20, minibatch="stochastic", ctl=list(batch_size=400))
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

### LOG-PCA
Sco <- CreateSeuratObject(counts=Y)
L <- median(colSums(Y))
YL <- log(sweep(Y,MARGIN=2,STATS=L^{-1}*colSums(Y),FUN="/") + 1)
my.pca <- irlba::prcomp_irlba(t(YL), n=20)
rownames(my.pca$x) <- colnames(Y)
Sco[["lpca"]] <- CreateDimReducObject(embeddings=my.pca$x,key="LPCA_")
Sco <- FindNeighbors(Sco,reduction = "lpca")
Sco <- FindClusters(Sco)
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
Sco <- FindClusters(Sco)
gbm <- adj.rand.index(phenoid, Sco$seurat_clusters)


#Umap
Sco <- RunUMAP(Sco, reduction = "gbm", dims=1:20)

p1 <- data.frame(x=Sco@reductions$umap@cell.embeddings[,1],
                 y=Sco@reductions$umap@cell.embeddings[,2],
                 color=phenoid) |>
  ggplot(aes(x=x,y=y,color=color)) +
  geom_point(size=0.5) + theme_bw() +
  xlab("UMAP1") + ylab("UMAP2") + labs(color="Cell type")

p2 <- data.frame(x=Sco@reductions$umap@cell.embeddings[,1],
                 y=Sco@reductions$umap@cell.embeddings[,2],
                 color=Sco$seurat_clusters) |>
  ggplot(aes(x=x,y=y,color=color)) +
  geom_point(size=0.5) + theme_bw() +
  xlab("UMAP1") + ylab("UMAP2") + labs(color="Cluster")

p.scgbm <- ggarrange(p1, p2, nrow=1) |>
  annotate_figure(top = text_grob("scGBM", face = "bold", size = 14))


set.seed(1)
fit <- glmpca(Y, L=20, minibatch="stochastic", ctl=list(batch_size=400))
Sco <- CreateSeuratObject(counts=Y)
Sco[["glmpca"]] <- CreateDimReducObject(embeddings=as.matrix(fit$res$factors),key="GLMPCA_")
Sco <- FindNeighbors(Sco,reduction = "glmpca")
Sco <- FindClusters(Sco)
glmpca <- adj.rand.index(phenoid, Sco$seurat_clusters)


Sco <- RunUMAP(Sco, dims=1:20, reduction="glmpca")

p1 <- data.frame(x=Sco@reductions$umap@cell.embeddings[,1],
                 y=Sco@reductions$umap@cell.embeddings[,2],
                 color=phenoid) |>
  ggplot(aes(x=x,y=y,color=color)) +
  geom_point(size=0.5) + theme_bw() +
  xlab("UMAP1") + ylab("UMAP2") + labs(color="Cell type")

p2 <- data.frame(x=Sco@reductions$umap@cell.embeddings[,1],
                 y=Sco@reductions$umap@cell.embeddings[,2],
                 color=Sco$seurat_clusters) |>
  ggplot(aes(x=x,y=y,color=color)) +
  geom_point(size=0.5) + theme_bw() +
  xlab("UMAP1") + ylab("UMAP2") + labs(color="Cluster")

p.glmpca <- ggarrange(p1, p2, nrow=1) |>
  annotate_figure(top = text_grob("GLM-PCA (SGD)", face = "bold", size = 14))


##Log+Scale+PCA
Sco <- CreateSeuratObject(counts=Y)
Sco <- NormalizeData(Sco)
Sco <- FindVariableFeatures(Sco,nfeatures=nrow(Y))
Sco <- ScaleData(Sco)
Sco <- RunPCA(Sco,npcs=20)
Sco <- FindNeighbors(Sco)
Sco <- FindClusters(Sco)
l2pca <- adj.rand.index(phenoid, Sco$seurat_clusters)

Sco <- RunUMAP(Sco, dims=1:20)

p1 <- data.frame(x=Sco@reductions$umap@cell.embeddings[,1],
                 y=Sco@reductions$umap@cell.embeddings[,2],
                 color=phenoid) |>
  ggplot(aes(x=x,y=y,color=color)) +
  geom_point(size=0.5) + theme_bw() +
  xlab("UMAP1") + ylab("UMAP2") + labs(color="Cell type")

p2 <- data.frame(x=Sco@reductions$umap@cell.embeddings[,1],
                 y=Sco@reductions$umap@cell.embeddings[,2],
                 color=Sco$seurat_clusters) |>
  ggplot(aes(x=x,y=y,color=color)) +
  geom_point(size=0.5) + theme_bw() +
  xlab("UMAP1") + ylab("UMAP2") + labs(color="Cluster")

p.lpca <- ggarrange(p1, p2, nrow=1) |>
  annotate_figure(top = text_grob("Log+Scale+PCA", face = "bold", size = 14))

##SCTRANSFORM
Sco <- CreateSeuratObject(counts=Y)
Sco <- SCTransform(Sco, variable.features.n = nrow(Y))
Sco <- RunPCA(Sco,npcs=20)
Sco <- FindNeighbors(Sco)
Sco <- FindClusters(Sco)
sct <- adj.rand.index(phenoid, Sco$seurat_clusters)

Sco <- RunUMAP(Sco, dims=1:20)

p1 <- data.frame(x=Sco@reductions$umap@cell.embeddings[,1],
                 y=Sco@reductions$umap@cell.embeddings[,2],
                 color=phenoid) |>
  ggplot(aes(x=x,y=y,color=color)) +
  geom_point(size=0.5) + theme_bw() +
  xlab("UMAP1") + ylab("UMAP2") + labs(color="Cell type")

p2 <- data.frame(x=Sco@reductions$umap@cell.embeddings[,1],
                 y=Sco@reductions$umap@cell.embeddings[,2],
                 color=Sco$seurat_clusters) |>
  ggplot(aes(x=x,y=y,color=color)) +
  geom_point(size=0.5) + theme_bw() +
  xlab("UMAP1") + ylab("UMAP2") + labs(color="Cluster")

p.sct <- ggarrange(p1, p2, nrow=1) |>
  annotate_figure(top = text_grob("SCT+PCA", face = "bold", size = 14))

###APR
Sco <- CreateSeuratObject(counts=Y)
Sco <- SCTransform(Sco, vst.flavor = "v1", method="offset",variable.features.n = nrow(Y))
Sco <- RunPCA(Sco,npcs=20)
Sco <- FindNeighbors(Sco)
Sco <- FindClusters(Sco)
#apr <- adj.rand.index(sce$phenoid, Sco$seurat_clusters)
apr <- adj.rand.index(phenoid, Sco$seurat_clusters)


Sco <- RunUMAP(Sco, dims=1:20)

p1 <- data.frame(x=Sco@reductions$umap@cell.embeddings[,1],
                 y=Sco@reductions$umap@cell.embeddings[,2],
                 color=phenoid) |>
  ggplot(aes(x=x,y=y,color=color)) +
  geom_point(size=0.5) + theme_bw() +
  xlab("UMAP1") + ylab("UMAP2") + labs(color="Cell type")

p2 <- data.frame(x=Sco@reductions$umap@cell.embeddings[,1],
                 y=Sco@reductions$umap@cell.embeddings[,2],
                 color=Sco$seurat_clusters) |>
  ggplot(aes(x=x,y=y,color=color)) +
  geom_point(size=0.5) + theme_bw() +
  xlab("UMAP1") + ylab("UMAP2") + labs(color="Cluster")

p.apr <- ggarrange(p1, p2, nrow=1) |>
  annotate_figure(top = text_grob("APR+PCA", face = "bold", size = 14))





##SCTRANSFORM
Sco <- CreateSeuratObject(counts=Y)
Sco <- SCTransform(Sco, variable.features.n = nrow(Y))
Sco <- RunPCA(Sco,npcs=20)
Sco <- FindNeighbors(Sco)
Sco <- FindClusters(Sco)
sct <- adj.rand.index(phenoid, Sco$seurat_clusters)

Sco <- RunUMAP(Sco, dims=1:20)

p1 <- data.frame(x=Sco@reductions$umap@cell.embeddings[,1],
                 y=Sco@reductions$umap@cell.embeddings[,2],
                 color=phenoid) |>
  ggplot(aes(x=x,y=y,color=color)) +
  geom_point(size=0.5) + theme_bw() +
  xlab("UMAP1") + ylab("UMAP2") + labs(color="Cell type")

p2 <- data.frame(x=Sco@reductions$umap@cell.embeddings[,1],
                 y=Sco@reductions$umap@cell.embeddings[,2],
                 color=Sco$seurat_clusters) |>
  ggplot(aes(x=x,y=y,color=color)) +
  geom_point(size=0.5) + theme_bw() +
  xlab("UMAP1") + ylab("UMAP2") + labs(color="Cluster")

p.sct <- ggarrange(p1, p2, nrow=1) |>
  annotate_figure(top = text_grob("SCT+PCA", face = "bold", size = 14))


### LOG-PCA
Sco <- CreateSeuratObject(counts=Y)
L <- median(colSums(Y))
YL <- log(sweep(Y,MARGIN=2,STATS=L^{-1}*colSums(Y),FUN="/") + 1)
my.pca <- irlba::prcomp_irlba(t(YL), n=20)
rownames(my.pca$x) <- colnames(Y)
Sco[["lpca"]] <- CreateDimReducObject(embeddings=my.pca$x,key="LPCA_")
Sco <- FindNeighbors(Sco,reduction = "lpca")
Sco <- FindClusters(Sco)
lpca_noscale <- adj.rand.index(phenoid, Sco$seurat_clusters)

Sco <- RunUMAP(Sco, dims=1:20, reduction="lpca")

p1 <- data.frame(x=Sco@reductions$umap@cell.embeddings[,1],
                 y=Sco@reductions$umap@cell.embeddings[,2],
                 color=phenoid) |>
  ggplot(aes(x=x,y=y,color=color)) +
  geom_point(size=0.5) + theme_bw() +
  xlab("UMAP1") + ylab("UMAP2") + labs(color="Cell type")

p2 <- data.frame(x=Sco@reductions$umap@cell.embeddings[,1],
                 y=Sco@reductions$umap@cell.embeddings[,2],
                 color=Sco$seurat_clusters) |>
  ggplot(aes(x=x,y=y,color=color)) +
  geom_point(size=0.5) + theme_bw() +
  xlab("UMAP1") + ylab("UMAP2") + labs(color="Cluster")

p.lpca_ns <- ggarrange(p1, p2, nrow=1) |>
  annotate_figure(top = text_grob("Log+PCA", face = "bold", size = 14))



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

results <- c(gbm,gbmproj, l2pca, sct,apr,glmpca,lpca_noscale)
names(results) <- c("scGBM-full", "scGBM-proj", "Log+Scale+PCA", "SCT", "APR","GLM-PCA (SGD)",
                    "Log+PCA")

saveRDS(results, file="../data/zhengmix8uneq.RDS")



##Final plot
library(ggpubr)

load(file="../data/scgbm_ggplot_obj.RData")
load(file="../data/glmpca_ggplot_obj.RData")
load(file="../data/scgbm_lpca_obj.RData")
load(file="../data/sct_ggplot_obj.RData")
load(file="../data/apr_ggplot_obj.RData")
load(file="../data/lpca_ns_ggplot_obj.RData")
load(file="../data/celltype_legend.RData")
load(file="../data/cluster_legend.RData")

p <- ggarrange(ggarrange(p.scgbm + guides(color="none"),
               p.glmpca,
               p.lpca,
               p.sct,
               p.apr,
               p.lpca_ns,
               nrow=3, ncol=2),
               ggarrange(celltype_legend,nrow=1),
               nrow=2,heights=c(10,1))



ggsave(p, filename="../plots/all_umap_subsampled.png", width=15, height=10, units="in")





##Resizing
load(file="../data/scgbm_ggplot_obj.RData")
load(file="../data/glmpca_ggplot_obj.RData")
load(file="../data/scgbm_lpca_obj.RData")
load(file="../data/sct_ggplot_obj.RData")
load(file="../data/apr_ggplot_obj.RData")
load(file="../data/lpca_ns_ggplot_obj.RData")





