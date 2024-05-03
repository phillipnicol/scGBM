

setwd("/Users/phillipnicol/Desktop/local_files/data/single-cell/Zheng_ERCC")
library(Seurat)
genes <- read.csv("genes.tsv")
barcodes <- read.csv("barcodes.tsv")

expr <- ReadMtx("matrix.mtx",cells="barcodes.tsv",
               features="genes.tsv")

Sco <-Seurat::CreateSeuratObject(counts=expr)

#Standard pipeline
Sco <- NormalizeData(Sco)
Sco <- FindVariableFeatures(Sco)
Sco <- ScaleData(Sco)
Sco <- RunPCA(Sco)
Sco <- RunUMAP(Sco,dims=1:10)
Sco <- FindNeighbors(Sco, dims = 1:10)
Sco <- FindClusters(Sco, resolution = 0.5)
DimPlot(Sco,group.by="seurat_clusters")
DimPlot(Sco,reduction="pca",group.by="seurat_clusters")


expr <- as.matrix(Sco@assays$RNA@counts)
expr <- expr[Sco@assays$RNA@var.features,]
mode(expr) <- "integer"


out <- gbm.nb(expr,M=1L)


A <- as.matrix(out$A)
B <- as.matrix(out$B)
c <- out$C[1,1]

for(i in 1:nrow(expr)) {
  for(j in 1:ncol(expr)) {
    Mu[i,j] <- out$A[i,1]+out$B[j,1]+out$C[1,1][[1]]
  }
}

Mu <- out$A


library(ggplot2)
V <- as.data.frame(out$V)
df <- data.frame(x=V[,1],y=V[,2],major=
                   as.character(clustering))

p <- ggplot(data=df,aes(x=x,y=y,color=major))
p <- p + geom_point()
p <- p + xlab("GBM 1")+ylab("GBM 2")
p <- p + theme_linedraw()
p

library(ggplot2)
V <- as.data.frame(out$V)
df <- data.frame(x0=V[,1],y0=V[,2],major=
                   "Monocyte",
                 a=1.96*out$se_V[,1],
                 b=1.96*out$se_V[,2])

p <- ggplot(data=df,aes(x0=x0,y0=y0,a=a,b=b,fill=major,
                        angle=0))
p <- p + geom_ellipse(alpha=0.5)
p <- p + xlab("GBM 1")+ylab("GBM 2")
p <- p + theme_linedraw()
p


u.max <- order(abs(out$U[,1]),decreasing=TRUE)[1:20]

df <- data.frame(u=out$U[u.max,1],
                 name=rownames(Y)[u.max])

p <- ggplot(data=df,aes(x=name,y=u))
p <- p + geom_point(color="blue")
p <- p + geom_hline(yintercept=0,color="red",lty="dashed")
p <- p + theme_classic()
p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p
