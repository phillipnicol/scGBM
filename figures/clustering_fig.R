
library(DuoClustering2018)
library(Seurat)
library(NewWave)
library(glmpca)
library(scGBM)
library(pdfCluster)

set.seed(1)
sce <- sce_full_Zhengmix8eq()

Sco <- as.Seurat(sce)
Sco <- FindVariableFeatures(Sco,nfeatures=2000)
Y <- Sco@assays$originalexp@counts[Sco@assays$originalexp@var.features,]
Y <- as.matrix(Y)

out <- gbm.sc(Y, max.iter=1000,tol=10^{-6}, M=20,infer.beta=TRUE)
#fit <- glmpca(Y,L=20,optimizer="avagrad",minibatch="stochastic",
#              ctl=list(batch_size=400,verbose=TRUE))

fit <- glmpca(Y,L=20)


V.glmpca <- as.matrix(fit$factors)

out <- get.se(out)

library(SingleCellExperiment)
my.sce <- SingleCellExperiment(list(counts=Y))
my.sce <- newWave(my.sce,K=20)

resolution <- seq(0.1,5,length.out=50)
df <- data.frame(gbm_standard = rep(0,length(resolution)),
                 gbm_coarse=rep(0,length(resolution)),
                 glmpca=rep(0,length(resolution)),
                 log = rep(0,length(resolution)),
                 sct= rep(0, length(resolution)),
                 nw = rep(0, length(resolution)))
df2 <- data.frame(gbm_standard = rep(0,length(resolution)),
                  gbm_coarse=rep(0,length(resolution)),
                  glmpca=rep(0,length(resolution)),
                  log = rep(0,length(resolution)),
                  sct= rep(0, length(resolution)),
                  nw = rep(0, length(resolution)))
combined <- sce$phenoid
combined[sce$phenoid %in% c("cd4.t.helper", "regulatory.t", "naive.t", "memory.t")] <- "t"
for(k in 1:length(resolution)) {
  Sco <- CreateSeuratObject(counts=Y)
  Sco[["gbm"]] <- CreateDimReducObject(embeddings=out$scores,key="GBM_")
  Sco <- FindNeighbors(Sco,reduction = "gbm")
  Sco <- FindClusters(Sco,resolution=resolution[k])

  df[k,1] <- adj.rand.index(Sco$seurat_clusters, sce$phenoid)
  df2[k,1] <- adj.rand.index(Sco$seurat_clusters, combined)
  cluster_fn <- function(V,Y) {
    Sco <- CreateSeuratObject(Y)
    colnames(V) <- 1:ncol(V)
    Sco[["gbm"]] <- CreateDimReducObject(embeddings=V,key="GBM_")
    Sco <- FindNeighbors(Sco,reduction = "gbm")
    Sco <- FindClusters(Sco,resolution=resolution[k])
    as.vector(Sco$seurat_clusters)
  }

  cci <- CCI(out,cluster.orig=Sco$seurat_clusters,reps=100,
             cluster.fn=cluster_fn,Y=Y)

  df[k,2] <- adj.rand.index(cci$coarse_cluster, sce$phenoid)
  df2[k,2] <- adj.rand.index(cci$coarse_cluster, combined)

  Sco <- CreateSeuratObject(counts=Y)
  Sco[["glmpca"]] <- CreateDimReducObject(embeddings=V.glmpca,key="GLMPCA_")
  Sco <- FindNeighbors(Sco,reduction = "glmpca")
  Sco <- FindClusters(Sco,resolution=resolution[k])
  df[k,3] <- adj.rand.index(Sco$seurat_clusters, sce$phenoid)
  df2[k,3] <- adj.rand.index(Sco$seurat_clusters, combined)

  Sco <- CreateSeuratObject(counts=Y)
  Sco <- NormalizeData(Sco)
  Sco <- FindVariableFeatures(Sco)
  Sco <- ScaleData(Sco)
  Sco <- RunPCA(Sco)
  Sco <- FindNeighbors(Sco)
  Sco <- FindClusters(Sco, resolution=resolution[k])
  df[k,4] <- adj.rand.index(Sco$seurat_clusters, sce$phenoid)
  df2[k,4] <- adj.rand.index(Sco$seurat_clusters, combined)

  Sco <- CreateSeuratObject(counts=Y)
  Sco <- SCTransform(Sco)
  Sco <- RunPCA(Sco)
  Sco <- FindNeighbors(Sco)
  Sco <- FindClusters(Sco, resolution=resolution[k])
  df[k,5] <- adj.rand.index(Sco$seurat_clusters, sce$phenoid)
  df2[k,5] <- adj.rand.index(Sco$seurat_clusters, combined)

  Sco <- CreateSeuratObject(counts=Y)
  Sco[["nw"]] <- CreateDimReducObject(embeddings=reducedDim(my.sce),key="NW_")
  Sco <- FindNeighbors(Sco,reduction = "nw")
  Sco <- FindClusters(Sco,resolution=resolution[k])
  df[k,6] <- adj.rand.index(Sco$seurat_clusters, sce$phenoid)
  df2[k,6] <- adj.rand.index(Sco$seurat_clusters, combined)

  print(df[k,])
  print(df2[k,])
}


library(ggplot2)
df$resolution <- resolution

df.gg <- reshape2::melt(df,id.vars="resolution")



p <- ggplot(data=df.gg,aes(x=resolution, y=value,
                           color=variable))
p <- p + geom_point()
p <- p + geom_line(linetype="dashed")
p <- p + xlab("Louvain resolution")
p <- p + ylab("ARI")
p <- p + theme_bw()


saveRDS(df,"results/8mix.RDS")
saveRDS(df2,"results/8mix2.RDS")

### Null DATA


#Clustering m
library(scGBM)
library(pheatmap)
library(ggpubr)
#Generate random noise
set.seed(1)
I <- 1000
J <- 5000

Y <- matrix(rpois(I*J,lambda=1),nrow=I,ncol=J)
rownames(Y) <- 1:I; colnames(Y) <- 1:J




out <- gbm.sc(Y, max.iter=1000, M=20,infer.beta=TRUE)
fit <- glmpca(Y,L=20,optimizer="avagrad",minibatch="stochastic",
              ctl=list(batch_size=400,verbose=TRUE))


V.glmpca <- as.matrix(fit$factors)

out <- get.se(out)

library(SingleCellExperiment)
my.sce <- SingleCellExperiment(list(counts=Y))
my.sce <- newWave(my.sce,K=20)


resolution <- seq(0.1,5,length.out=50)
df <- data.frame(gbm_standard = rep(0,length(resolution)),
                 gbm_coarse=rep(0,length(resolution)),
                 glmpca=rep(0,length(resolution)),
                 log = rep(0,length(resolution)),
                 sct= rep(0, length(resolution)),
                 nw = rep(0, length(resolution)))
for(k in 1:length(resolution)) {
  Sco <- CreateSeuratObject(counts=Y)
  Sco[["gbm"]] <- CreateDimReducObject(embeddings=out$scores,key="GBM_")
  Sco <- FindNeighbors(Sco,reduction = "gbm")
  Sco <- FindClusters(Sco,resolution=resolution[k])

  df[k,1] <- length(table(Sco$seurat_clusters))
  cluster_fn <- function(V,Y) {
    Sco <- CreateSeuratObject(Y)
    colnames(V) <- 1:ncol(V)
    Sco[["gbm"]] <- CreateDimReducObject(embeddings=V,key="GBM_")
    Sco <- FindNeighbors(Sco,reduction = "gbm")
    Sco <- FindClusters(Sco,resolution=resolution[k])
    as.vector(Sco$seurat_clusters)
  }

  cci <- CCI(out,cluster.orig=Sco$seurat_clusters,reps=100,
             cluster.fn=cluster_fn,Y=Y)

  df[k,2] <-  length(table(cci$coarse_cluster))

  Sco <- CreateSeuratObject(counts=Y)
  Sco[["glmpca"]] <- CreateDimReducObject(embeddings=V.glmpca,key="GLMPCA_")
  Sco <- FindNeighbors(Sco,reduction = "glmpca")
  Sco <- FindClusters(Sco,resolution=resolution[k])
  df[k,3] <- length(table(Sco$seurat_clusters))

  Sco <- CreateSeuratObject(counts=Y)
  Sco <- NormalizeData(Sco)
  Sco <- FindVariableFeatures(Sco)
  Sco <- ScaleData(Sco)
  Sco <- RunPCA(Sco)
  Sco <- FindNeighbors(Sco)
  Sco <- FindClusters(Sco, resolution=resolution[k])
  df[k,4] <-  length(table(Sco$seurat_clusters))

  Sco <- CreateSeuratObject(counts=Y)
  Sco <- SCTransform(Sco)
  Sco <- RunPCA(Sco)
  Sco <- FindNeighbors(Sco)
  Sco <- FindClusters(Sco, resolution=resolution[k])
  df[k,5] <-  length(table(Sco$seurat_clusters))

  Sco <- CreateSeuratObject(counts=Y)
  Sco[["nw"]] <- CreateDimReducObject(embeddings=reducedDim(my.sce),key="NW_")
  Sco <- FindNeighbors(Sco,reduction = "nw")
  Sco <- FindClusters(Sco,resolution=resolution[k])
  df[k,6] <- length(table(Sco$seurat_clusters))

  print(df[k,])
}



df$resolution <- resolution

df.gg <- reshape2::melt(df,id.vars="resolution")



p <- ggplot(data=df.gg,aes(x=resolution, y=value,
                           color=variable))
p <- p + geom_point()
p <- p + geom_line(linetype="dashed")
p <- p + xlab("Louvain resolution")
p <- p + ylab("# of clusters")
p <- p + theme_bw()

saveRDS(df,"results/null.RDS")










## Directory for data
data <- 0 ## TO DO PUT YOUR DIRECTORY HERE

#Clustering m
library(scGBM)
library(pheatmap)
library(ggpubr)

#Generate random noise
set.seed(1)
I <- 1000
J <- 5000

pt.size <- 0.5

Y <- matrix(rpois(I*J,lambda=1),nrow=I,ncol=J)
rownames(Y) <- 1:I; colnames(Y) <- 1:J

library(Seurat)

#Sco <- CreateSeuratObject(counts=Y)
#Sco <- NormalizeData(Sco)
#Sco <- FindVariableFeatures(Sco)
#Sco <- ScaleData(Sco)
#Sco <- RunPCA(Sco)

#Sco <- FindNeighbors(Sco)
#Sco <- FindClusters(Sco)


# Fit scGBM

out <- gbm.sc(Y,M=20)

out <- get.se(out)


colnames(Y) <- 1:5000; rownames(Y) <- 1:1000
Sco <- CreateSeuratObject(counts=Y)
V <- out$scores; colnames(out$scores) <- 1:20
Sco[["gbm"]] <- CreateDimReducObject(embeddings=out$V,key="GBM_")
Sco <- FindNeighbors(Sco,reduction = "gbm")
Sco <- FindClusters(Sco)

cluster_fn <- function(V,Y) {
  Sco <- CreateSeuratObject(Y)
  colnames(V) <- 1:ncol(V)
  Sco[["gbm"]] <- CreateDimReducObject(embeddings=V,key="GBM_")
  Sco <- FindNeighbors(Sco,reduction = "gbm")
  Sco <- FindClusters(Sco)
  as.vector(Sco$seurat_clusters)
}

cci.1 <- CCI(out,cluster.orig=Sco$seurat_clusters,
             cluster.fn = cluster_fn, Y=Y, null.dist = TRUE)

save(out, file=paste0(data, "Miller/Clustering_fig/null_scGBM.RData"))
saveRDS(Sco$seurat_clusters, paste0(data, "Miller/Clustering_fig/null_clusterorig.RData"))
save(cci.1, file=paste0(data, "Miller/Clustering_fig/null_cci.RData"))



### ZHENG mix data

library(DuoClustering2018)
library(Seurat)
library(NewWave)
library(glmpca)
library(scGBM)
library(pdfCluster)

set.seed(1)
sce <- sce_full_Zhengmix8eq()

Sco <- as.Seurat(sce)
Sco <- FindVariableFeatures(Sco,nfeatures=2000)
Y <- Sco@assays$originalexp@counts[Sco@assays$originalexp@var.features,]
Y <- as.matrix(Y)

out <- gbm.sc(Y, max.iter=1000,tol=10^{-6}, M=20,infer.beta=TRUE)

out <- get.se(out, EPS=0.001)


cci <- CCI(out,cluster.orig=sce$phenoid,
           centers=8,
           null.dist = TRUE)

save(out, file=paste0(data.dir ,"8mix_scGBM.Rdata"))

save(cci, file=paste0(data.dir, "8mix_cci.RData"))



### NULL DATA

## Directory for data
data <- "/Volumes/pnicol_data/"

load(paste0(data, "null_scGBM.RData"))
load(paste0(data, "null_cci.RData"))
cluster.orig <- readRDS(paste0(data, "null_clusterorig.RData"))

library(ggplot2); library(ggforce)
library(ggpubr)
library(scGBM)
pa1 <- plot_gbm(out,se=TRUE,return.gg=TRUE,cluster=cluster.orig) + theme_bw()
pa1 <- pa1 + ggtitle("Null data") + xlab("")

pb1 <- cci.1$cci_diagonal + ggtitle("Distribution of Cluster Confidence Index (CCI)") + ylim(0,1)
pb1 <- pb1 + xlab("")

H.table <- cci.1$H.table
rownames(H.table) <- unique(cluster.orig)
colnames(H.table) <- unique(cluster.orig)
pc1 <- pheatmap::pheatmap(H.table,color=colorRampPalette(c("white","red"))(100),
                          breaks=seq(0,1,by=0.01),
                          cluster_cols=FALSE,
                          rownames=TRUE,
                          colnames=TRUE,
                          cluster_rows = FALSE,
                          legend=FALSE,
                          treeheight_row = 0,
                          treeheight_col = 0)[[4]]


nullcoarse <- readRDS(paste0(data, "null.RDS"))

df.gg <- reshape2::melt(nullcoarse,id.vars="resolution")



p <- ggplot(data=df.gg,aes(x=resolution, y=value,
                           color=variable))
p <- p + geom_point()
p <- p + geom_line(linetype="dashed")
p <- p + xlab("Louvain resolution")
p <- p + ylab("# of clusters")
p <- p + scale_color_hue(labels=c("scGBM",
                                  "scGBM+CCI",
                                  "GLM-PCA",
                                  "Log+PCA",
                                  "SCT+PCA",
                                  "NewWave"))
p <- p + labs(color="Method")
p <- p + theme_bw()
pd <- p + xlim(0,4)
ggsave(pd, filename=paste0(data,"null_coarse.pdf"),units="in",width=8, height=4)


### Zheng Data

## Directory for data
data <- "/Volumes/pnicol_data/"

library(DuoClustering2018)
sce <- sce_full_Zhengmix8eq()


load(paste0(data, "8mix_scGBM.RData"))
load(paste0(data, "8mix_cci.RData"))

library(ggplot2); library(ggforce)
shortened.cluster <- sce$phenoid
shortened.cluster[sce$phenoid == "cd14.monocytes"] <- "M"
shortened.cluster[sce$phenoid == "b.cells"] <- "B"
shortened.cluster[sce$phenoid == "naive.cytotoxic"] <- "NC"
shortened.cluster[sce$phenoid == "regulatory.t"] <- "rT"
shortened.cluster[sce$phenoid == "cd4.t.helper"] <- "Th"
shortened.cluster[sce$phenoid == "cd56.nk"] <- "NK"
shortened.cluster[sce$phenoid == "memory.t"] <- "mT"
shortened.cluster[sce$phenoid == "naive.t"] <- "nT"
pa2 <- plot_gbm(out,se=TRUE,return.gg=TRUE,cluster=shortened.cluster) + theme_bw()
pa2 <- pa2 + ggtitle("10X immune cells")


pb2 <- cci$cci_diagonal + ggtitle("")+ ylim(0,1)
pb2 <- pb2 + scale_x_discrete(labels=c("B", "NC", "M", "rT", "Th", "NK", "mT", "nT"))
pb2 <- pb2 + ylab("") + xlab("")


H.table <- cci$H.table
rownames(H.table) <- c("B", "NC", "M", "rT", "Th", "NK", "mT", "nT")
colnames(H.table) <- c("B", "NC", "M", "rT", "Th", "NK", "mT", "nT")
pc2 <- pheatmap::pheatmap(H.table,legend=TRUE, color=colorRampPalette(c("white","red"))(100),
                          breaks=seq(0,1,by=0.01),
                          cluster_cols=TRUE,
                          rownames=TRUE,
                          colnames=TRUE,
                          cluster_rows = TRUE,
                          treeheight_row = 0,
                          treeheight_col = 0)[[4]]


emixcoarse <- readRDS(paste0(data, "8mix.RDS"))
emix2 <- readRDS(paste0(data, "8mix2.RDS"))

df.gg <- reshape2::melt(emixcoarse,id.vars="resolution")




p <- ggplot(data=rev(df.gg),aes(x=resolution, y=value,
                                color=variable))
p <- p + geom_point()
p <- p + geom_line(linetype="dashed")
p <- p + xlab("Louvain resolution")
p <- p + ylab("ARI")
p <- p + scale_color_hue(labels=c("scGBM",
                                  "scGBM+CCI",
                                  "GLM-PCA",
                                  "Log+PCA",
                                  "SCT+PCA",
                                  "NewWave"))
p <- p + labs(color="Method")
p <- p + theme_bw()
pd1 <- p + xlim(0,4)


emix2$resolution <- emixcoarse$resolution
df.gg <- reshape2::melt(emix2,id.vars="resolution")



p <- ggplot(data=df.gg,aes(x=resolution, y=value,
                           color=variable))
p <- p + geom_point()
p <- p + geom_line(linetype="dashed")
p <- p + xlab("Louvain resolution")
p <- p + ylab("")
p <- p + scale_color_hue(labels=c("scGBM",
                                  "scGBM+CCI",
                                  "GLM-PCA",
                                  "Log+PCA",
                                  "SCT+PCA",
                                  "NewWave"))
p <- p + labs(color="Method")
pd2 <- p + theme_bw() + xlim(0,4)

pbot <- ggarrange(pd1, pd2, common.legend = TRUE,
                  legend = "bottom")

p <- ggarrange(ggarrange(pa1,pa2,pb1,pb2,pc1,pc2,ncol=2,nrow=3,
                         labels=c("a", "","b","","c","")),
               ggarrange(pd1, pd2, common.legend = TRUE,
                         legend = "bottom",nrow=1,
                         labels=c("d","e")),nrow=2,ncol=1,
               heights=c(7,3))
p


