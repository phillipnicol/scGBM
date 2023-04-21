


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

out <- gbm.sc(Y,M=10)

library(ggplot2)

colnames(Y) <- 1:5000; rownames(Y) <- 1:1000
Sco <- CreateSeuratObject(counts=Y)
V <- out$V; colnames(out$V) <- 1:10
Sco[["gbm"]] <- CreateDimReducObject(embeddings=out$V,key="GBM_")
Sco <- FindNeighbors(Sco,reduction = "gbm")
Sco <- FindClusters(Sco)

df <- data.frame(V1=out$V[,1],V2=out$V[,2],
                 color=Sco$seurat_clusters)
p <- ggplot(data=df,aes(x=V1,y=V2,color=color))
p <- p + geom_point(size=pt.size)
p <- p + xlab("GBM-1")+ylab("GBM-2") +  theme_bw()
p <- p + labs(color="Cluster")
p_A <- p

# Uncertainty orbs
library(ggforce)
out <- get.se(out)
p <- plot.gbm(out,cluster=Sco$seurat_clusters,
              se=TRUE,return.gg=TRUE)
p <- p + xlab("GBM-1")+ylab("GBM-2")+theme_bw()+theme(legend.position="none")
p_B <- p


cluster_fn <- function(V,Y) {
  Sco <- CreateSeuratObject(Y)
  colnames(V) <- 1:ncol(V)
  Sco[["gbm"]] <- CreateDimReducObject(embeddings=V,key="GBM_")
  Sco <- FindNeighbors(Sco,reduction = "gbm")
  Sco <- FindClusters(Sco)
  as.vector(Sco$seurat_clusters)
}

cci.1 <- CCI(out,cluster.orig=Sco$seurat_clusters,
             cluster.fn = cluster_fn, Y=Y)
p_C <- cci.1$cci.plot
