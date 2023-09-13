data.dir <- 0 ##PUT YOUR PATH HERE

library(Seurat)
library(glmpca)
library(scGBM)
set.seed(1)
points.size <- 0.5
I <- 1000
J <- 1000
Mu <- matrix(1,nrow=I,ncol=J)
Mu[1,333:666] <- 10
Mu[2,667:1000] <- 50
Y <- matrix(rpois(n=I*J,lambda=as.vector(Mu)),nrow=I,ncol=J)

true_cluster <- rep(1,J)
true_cluster[1:332] <- "C"
true_cluster[333:666] <- "A"
true_cluster[667:1000] <- "B"
true_cluster <- as.character(true_cluster)
colnames(Y) <- 1:J
rownames(Y) <- 1:I
Sco <- CreateSeuratObject(counts=Y)
Sco$group <- true_cluster

Sco <- NormalizeData(Sco) #Verified log(1+CPT, base=e) transform
Sco <- FindVariableFeatures(Sco)
Sco <- ScaleData(Sco)
Sco <- RunPCA(Sco)
lpca <- Sco@reductions$pca@cell.embeddings
df <- data.frame(x=lpca[,1],y=lpca[,2],color=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=points.size)
p <- p + theme_bw()+xlab("PC1")+ylab("PC2")+labs(color="Cell Type")
p <- p + ggtitle("log+scale+PCA")
p_log2 <- p
Sco <- RunUMAP(Sco,dims=1:20,reduction="pca")
lpca.umap <- Sco@reductions$umap@cell.embeddings
Sco <- RunTSNE(Sco)
lpca.tsne <- Sco@reductions$tsne@cell.embeddings
df <- data.frame(x=lpca.umap[,1],y=lpca.umap[,2],color=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=points.size)
p <- p + theme_bw()+xlab("UMAP 1")+ylab("UMAP 2")+labs(color="Cell Type")
p <- p + ggtitle("log+scale+PCA+UMAP")
p_log2umap <- p


U <- Sco@reductions$pca@feature.loadings[,1]
df <- data.frame(x=1:I,y=U,color=ifelse(1:I<3,"red","grey"))
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=points.size)
p <- p + scale_color_manual(values=c("grey","red"))
p <- p + theme_bw()+xlab("Gene")+ylab("Factor 1 weight")+guides(color="none")
p <- p + ggtitle("log+scale+PCA")
p_log2U <- p

Sco <- SCTransform(Sco)
Sco <- RunPCA(Sco)
sct <- Sco@reductions$pca@cell.embeddings
df <- data.frame(x=sct[,1],y=sct[,2],color=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=points.size)
p <- p + theme_bw()+xlab("PC1")+ylab("PC2")+guides(color="none")
p <- p + ggtitle("SCT+PCA")
p_sct <- p
Sco <- RunUMAP(Sco,dims=1:20)
sct.umap <- Sco@reductions$umap@cell.embeddings
Sco <- RunTSNE(Sco)
sct.tsne <- Sco@reductions$tsne@cell.embeddings

df <- data.frame(x=sct.umap[,1],y=sct.umap[,2],color=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=points.size)
p <- p + theme_bw()+xlab("UMAP 1")+ylab("UMAP 2")+labs(color="Cell Type")
p <- p + ggtitle("SCT+PCA+UMAP")+guides(color="none")
p_sctumap <- p

U <- Sco@reductions$pca@feature.loadings[,1]
df <- data.frame(x=1:I,y=U,color=ifelse(1:I<3,"red","grey"))
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=points.size)
p <- p + scale_color_manual(values=c("grey","red"))
p <- p + theme_bw()+xlab("Gene")+ylab("Factor 1 weight")+guides(color="none")
p <- p + ggtitle("SCT+PCA")
p_sctU <- p

p <- ggarrange(p_sct, p_log2, nrow=1,ncol=2)
ggsave(filename=paste0(data.dir, "sct_pca_scores_mixture_cells.png"),
       units="in", width=11,height=4)

out <- gbm.sc(Y,M=20,tol=10^{-10}, infer.beta=TRUE)

Sco <- CreateSeuratObject(counts=Y)
Sco[["gbm"]] <- CreateDimReducObject(embeddings=out$scores,key="GBM_")
Sco <- RunUMAP(Sco,reduction="gbm",dims=1:20)
df <- data.frame(x=Sco@reductions$umap@cell.embeddings[,1],
                 y=Sco@reductions$umap@cell.embeddings[,2],color=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=points.size)
p <- p + theme_bw() + guides(color="none")
p <- p+xlab("GBM1")+ylab("GBM2")+ggtitle("scGBM")
p_gbmUMAP <- p


df <- data.frame(x=out$V[,1],y=out$V[,2],color=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=points.size)
p <- p + theme_bw() + guides(color="none")
p <- p+xlab("GBM1")+ylab("GBM2")+ggtitle("scGBM")
p_gbm <- p

#U <- Sco@reductions$pca@feature.loadings[,1]
df <- data.frame(x=1:I,y=out$U[,1],color=ifelse(1:I<3,"red","grey"))
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=points.size)
p <- p + scale_color_manual(values=c("grey","red"))
p <- p + theme_bw()+xlab("Gene")+ylab("Factor 1 weight")+guides(color="none")
p <- p + ggtitle("scGBM")
p_gbmU <- p

df <- data.frame(y=Y[1,],x=true_cluster,fill=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,fill=fill))+geom_boxplot()
p <- p + theme_bw() + guides(fill="none")
p <- p + xlab("Cluster") + ylab("Counts")
p_g1 <- p

df <- data.frame(y=Y[2,],x=true_cluster,fill=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,fill=fill))+geom_boxplot()
p <- p + theme_bw()+ guides(fill="none")
p <- p + xlab("Cluster") + ylab("Counts")
p_g2 <- p

df <- data.frame(y=Y[3,],x=true_cluster,fill=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,fill=fill))+geom_boxplot()
p <- p + theme_bw()+guides(fill="none")
p <- p + xlab("Cluster") + ylab("Counts")
p_g3 <- p

df <- data.frame(x=lpca.tsne[,1],y=lpca.tsne[,2],color=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=points.size)
p <- p + theme_bw()+xlab("T1")+ylab("T2")+labs(color="Cell Type")
p <- p + ggtitle("log+scale+PCA+tSNE")
p_logtsne <- p

df <- data.frame(x=sct.tsne[,1],y=sct.tsne[,2],color=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=points.size)
p <- p + theme_bw()+xlab("T1")+ylab("T2")+labs(color="Cell Type")
p <- p + ggtitle("SCT+PCA+tSNE")
p_scttsne <- p


p <- ggarrange(p_logtsne,p_scttsne,
               nrow=1,ncol=2, common.legend=TRUE,
               legend = "bottom")
ggsave(p, filename=paste0(data.dir,"tnse_plot.png"),
       width=6.5, height=3, units="in")

out50 <- gbm.sc(Y, M=50, tol=10^{-10}, infer.beta=TRUE)
cor(out$scores[,1], out50$scores[,1])^2
cor(out$scores[,2], out50$scores[,2])^2
#RUN GLM-PCA 6 times
pt <- list()
library(glmpca)
for(t in 1:6) {
  print(t)
  fit <- glmpca(Y, L=50)

  U <- fit$loadings
  V <- fit$factors

  df <- data.frame(x=V[,1],y=V[,2],color=true_cluster)
  p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=points.size)
  p <- p + theme_bw() + guides(color="none")
  p <- p+xlab("V1")+ylab("V2")
  p_glmpca <- p
  pt[[t]] <- p
}


ggsave(pt[[1]], filename=paste0(data.dir,"sim_comparison_glmpca.png"), units="in", width=11,height=4)





### Download monocyte data from 10X
library(Seurat)
genes <- read.csv("genes.tsv")
barcodes <- read.csv("barcodes.tsv")

expr <- ReadMtx("matrix.mtx",cells="barcodes.tsv",
                features="genes.tsv")

Sco <-Seurat::CreateSeuratObject(counts=expr)

Sco <- NormalizeData(Sco)
Sco <- FindVariableFeatures(Sco)
Sco <- ScaleData(Sco)
Sco <- RunPCA(Sco)

Y <- Sco@assays$RNA@counts[Sco@assays$RNA@var.features,]
Y <- as.matrix(Y)
Y <- Y[rowSums(Y) >= 10,]

out <- gbm.sc(Y,M=20,max.iter=1000)
library(ashr)
out <- scGBM::get.se(out)
U.correct <- get_pm(ash(out$U[,1],out$se_U[,1]))
gene.means <- rowMeans(Y)
df <- data.frame(x=log(gene.means),y=out$U[,1],color="Raw")
df2 <- data.frame(x=log(gene.means),y=U.correct,color="Stabilized")
df <- rbind(df,df2)
p <- ggplot(data=df,aes(x=x,y=y,color=color))
p <- p + geom_point(size=points.size)
p <- p + scale_color_manual(values=c("grey","red"),
                            labels=c("Raw","Stabilized"))
p <- p+xlab("Log mean count")+ylab("Factor 1 weight")+labs(color="")
p <- p + theme_bw()
p_mono <- p
p_mono <- p_mono + ggtitle("Monocytes")
