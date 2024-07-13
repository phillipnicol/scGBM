
library(tidyverse)
set.seed(1)
library(Seurat)
library(viridis)
pt.size <- 0.5

set.seed(1)
I <- 1000
J <- 1000
#baseline.means <- rexp(n=I, rate=1)
baseline.means <- rep(1, I)
Mu <- matrix(baseline.means,nrow=I,ncol=J)
Mu[1,] <- 100
Mu[1,1:10] <- 1
Mu[1,11:20] <- 10
Mu[2,] <- 1
Mu[2,667:1000] <- 50
Y <- matrix(rpois(n=I*J,lambda=as.vector(Mu)),nrow=I,ncol=J)


true_cluster <- rep(1,J)
true_cluster[1:10] <- "A"
true_cluster[11:20] <- "B"
true_cluster[21:666] <- "C"
true_cluster[667:1000] <- "D"
true_cluster <- as.character(true_cluster)

df <- data.frame(y=Y[1,],x=true_cluster,fill=true_cluster)
pg1 <- ggplot(data=df,aes(x=x,y=y,fill=fill))+geom_boxplot()+
  theme_bw() + guides(fill="none") +
  geom_jitter(alpha=0.5, size=0.5,width=0.3) +
  xlab("") + ylab("Counts") + ggtitle("Gene 1") +
  scale_fill_manual(values = c("A" = "#FF0000", # Bright red
                               "B" = "#0000FF", # Bright blue
                               "C" = "#CCCCCC", # Light grey
                               "D" = "#999999"))+  # Darker grey
  scale_y_sqrt()

df <- data.frame(y=Y[2,],x=true_cluster,fill=true_cluster)
pg2 <- ggplot(data=df,aes(x=x,y=y,fill=fill))+geom_boxplot()+
  theme_bw() + guides(fill="none") +
  geom_jitter(alpha=0.5, size=0.5,width=0.3) +
  xlab("") + ylab("Counts") + ggtitle("Gene 2") +
  scale_fill_manual(values = c("A" = "#FF0000", # Bright red
                                "B" = "#0000FF", # Bright blue
                                "C" = "#CCCCCC", # Light grey
                                "D" = "#999999"))+  # Darker grey
  scale_y_sqrt()

pg3 <- Y[3:1000,] |> rbind(true_cluster) |>
  t() |>
  as.data.frame() |>
  pivot_longer(cols=-c(999)) |>
  mutate(value=as.numeric(value)) |>
  sample_n(10^4) |>
  ggplot(aes(x=true_cluster,y=value,fill=true_cluster))+geom_boxplot()+
  theme_bw() + guides(fill="none") +
  geom_jitter(alpha=0.5, size=0.5,width=0.3) +
  xlab("") + ylab("Counts") + ggtitle("Genes 3-1000 (random noise)") +
  scale_fill_manual(values = c("A" = "#FF0000", # Bright red
                               "B" = "#0000FF", # Bright blue
                               "C" = "#CCCCCC", # Light grey
                               "D" = "#999999"))+  # Darker grey
  scale_y_sqrt()


### LOG + SCALE + PCA (SEURAT)
colnames(Y) <- 1:J
rownames(Y) <- 1:I
Sco <- CreateSeuratObject(counts=Y)
Sco <- NormalizeData(Sco)
Sco <- FindVariableFeatures(Sco)
Sco <- ScaleData(Sco)
Sco$group <- true_cluster
Sco <- RunPCA(Sco)
lpca <- Sco@reductions$pca@cell.embeddings
lpca.scale.umap <- umap::umap(lpca[,1:20])$layout
df <- data.frame(x=lpca[,1],y=lpca[,2],color=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=pt.size)
p <- p + theme_bw()+xlab("")+ylab("")+guides(color="none") +
  ggtitle("Log+Scale+PCA") + theme(plot.title = element_text(size = 10)) +
  scale_color_manual(values = c("A" = "#FF0000", # Bright red
                                "B" = "#0000FF", # Bright blue
                                "C" = "#CCCCCC", # Light grey
                                "D" = "#999999"))  # Darker grey
p_lpcascale <- p

df <- data.frame(x=lpca.scale.umap[,1],
                 y=lpca.scale.umap[,2],color=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=pt.size)
p <- p + theme_bw()+xlab("")+ylab("")+guides(color="none") +
  ggtitle("LOG+SCALE+PCA+UMAP") +
  scale_color_manual(values = c("A" = "#FF0000", # Bright red
                                "B" = "#0000FF", # Bright blue
                                "C" = "#CCCCCC", # Light grey
                                "D" = "#999999"))  # Darker grey
p_lpcascale_umap <- p




### SCTRANSFORM (SEURAT)
colnames(Y) <- 1:J
rownames(Y) <- 1:I
Sco <- CreateSeuratObject(counts=Y)
Sco$group <- true_cluster
Sco <- SCTransform(Sco)
Sco <- RunPCA(Sco)
sct <- Sco@reductions$pca@cell.embeddings
sct.umap <- umap::umap(sct[,1:20])$layout
df <- data.frame(x=sct[,1],y=sct[,2],color=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=pt.size)
p <- p + theme_bw()+xlab("")+ylab("")+guides(color="none") +
  ggtitle("SCT+PCA") +
  scale_color_manual(values = c("A" = "#FF0000", # Bright red
                                "B" = "#0000FF", # Bright blue
                                "C" = "#CCCCCC", # Light grey
                                "D" = "#999999"))  # Darker grey
p_sct <- p

df <- data.frame(x=sct.umap[,1],y=sct.umap[,2],color=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=pt.size)
p <- p + theme_bw()+xlab("")+ylab("")+guides(color="none") +
  ggtitle("SCT+PCA+UMAP") +
  scale_color_manual(values = c("A" = "#FF0000", # Bright red
                                "B" = "#0000FF", # Bright blue
                                "C" = "#CCCCCC", # Light grey
                                "D" = "#999999"))  # Darker grey
p_sct_umap <- p

g1 <- Sco@assays$SCT@scale.data[1,]
g2 <- Sco@assays$SCT@scale.data[2,]
nullg <- as.vector(Sco@assays$SCT@scale.data[3:1000,])
df <- data.frame(x=rep(c("Gene 1", "Gene 2", "Null Genes"), times=c(1000, 1000, 998000)),
                 y=c(g1,g2,nullg))
p <- ggplot(data=df,aes(x=x,y=abs(y),fill=x)) +
  geom_boxplot()
p




### ANALYTIC PEARSON RESIDUALS
apr <- sctransform::vst(Y, method="offset")
pca.apr <- prcomp(t(apr$y))
apr.umap <- umap::umap(pca.apr$x[,1:20])$layout
df <- data.frame(x=pca.apr$x[,1],y=pca.apr$x[,2],color=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=pt.size)
p <- p + theme_bw()+xlab("")+ylab("")+guides(color="none")+
  ggtitle("APR+PCA") +
  scale_color_manual(values = c("A" = "#FF0000", # Bright red
                                "B" = "#0000FF", # Bright blue
                                "C" = "#CCCCCC", # Light grey
                                "D" = "#999999"))  # Darker grey
p_apr <- p

df <- data.frame(x=apr.umap[,1],y=apr.umap[,2],color=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=pt.size)
p <- p + theme_bw()+xlab("")+ylab("")+guides(color="none") +
  ggtitle("APR+PCA+UMAP") +
  scale_color_manual(values = c("A" = "#FF0000", # Bright red
                                "B" = "#0000FF", # Bright blue
                                "C" = "#CCCCCC", # Light grey
                                "D" = "#999999"))  # Darker grey
p_apr_umap <- p

gene.var.apr <- apply(apr$y, 1, var)[-2]
df <- data.frame(x=rep("APR",998),
                 y=gene.var.apr[-1])

p <- ggplot(data=df,aes(x=x,y=y)) + geom_boxplot(outlier.shape = NA,fill="lightblue")+
  geom_jitter(width=0.2,size=0.5) +
  geom_point(x="APR",y=gene.var.apr[1], color="red", size=2, shape=8) +
  ylim(0.8,1.35) +
  xlab("") + ylab("Variance") + theme_bw()

### LOG +  PCA
colnames(Y) <- 1:J
rownames(Y) <- 1:I
L <- median(colSums(Y))
YL <- log(sweep(Y,MARGIN=2,STATS=L^{-1}*colSums(Y),FUN="/") + 1)
my.pca <- prcomp(t(YL))
lpca <- my.pca$x
lpca.umap <- umap::umap(lpca[,1:20])$layout
df <- data.frame(x=lpca[,1],y=lpca[,2],color=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=pt.size)
p <- p + theme_bw()+xlab("")+ylab("")+guides(color="none") +
  ggtitle("Log+PCA") +
  scale_color_manual(values = c("A" = "#FF0000", # Bright red
                               "B" = "#0000FF", # Bright blue
                               "C" = "#CCCCCC", # Light grey
                               "D" = "#999999"))  # Darker grey
p_lpca <- p

df <- data.frame(x=lpca.umap[,1],y=lpca.umap[,2],color=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=pt.size)
p <- p + theme_bw()+xlab("")+ylab("")+guides(color="none") +
  ggtitle("APR+PCA+UMAP") +
  scale_color_manual(values = c("A" = "#FF0000", # Bright red
                                "B" = "#0000FF", # Bright blue
                                "C" = "#CCCCCC", # Light grey
                                "D" = "#999999"))  # Darker grey
p_lpca_umap <- p

gene.var.log <- apply(YL, 1, var)[-2]
df <- data.frame(x=rep("LOG",998),
                 y=gene.var.log[-1])

p <- ggplot(data=df,aes(x=x,y=y)) + geom_boxplot(outlier.shape = NA,fill="lightblue")+
  geom_jitter(width=0.2,size=0.5) +
  geom_point(x="LOG",y=gene.var.log[1], color="red", size=2, shape=8) +
  xlab("") + ylab("Variance") + theme_bw()
p <- p_log_var

library(ggpubr)
p.single.full <- ggarrange(ggarrange(pg1, pg2, pg3, nrow=1),
          ggarrange(p_lpca, p_lpcascale,
                    p_sct, p_apr, nrow=1), nrow=2,
        heights=c(1,1))
ggsave(p.single.full,
       filename="../plots/single_marker_pca.png",
       width=11.9, height=6.32, units="in")


p.umap <- ggarrange(p_lpca_umap, p_lpcascale_umap,
                    p_sct_umap, p_apr_umap)


###scGBM
set.seed(42)
out <- gbm.sc(Y,M=20,sigma=10, infer.beta=TRUE)


df <- data.frame(x=out$scores[,1], y=out$scores[,2],color=true_cluster)

p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=pt.size)
p <- p + theme_bw()+xlab("GBM1")+ylab("GBM1")+guides(color="none") +
  scale_color_manual(values = c("A" = "#FF0000", # Bright red
                              "B" = "#0000FF", # Bright blue
                              "C" = "#CCCCCC", # Light grey
                              "D" = "#999999")) +
  ggtitle("Simulated data with four clusters")

ggsave(p, filename="../plots/single_marker.png")


df <- data.frame(x=my.umap$layout[,1],y=my.umap$layout[,2],color=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=pt.size)
p <- p + theme_bw()+xlab("")+ylab("")+guides(color="none") +
  ggtitle("Log+PCA")
p_umap <- p

