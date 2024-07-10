
library(tidyverse)
set.seed(1)
library(Seurat)
library(viridis)
pt.size <- 0.5
L <- 10^4

set.seed(1)
I <- 1000
J <- 1000
baseline.means <- rexp(n=I, rate=1)
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
df <- data.frame(x=lpca[,1],y=lpca[,2],color=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=pt.size)
p <- p + theme_bw()+xlab("")+ylab("")+guides(color="none")+
  ggtitle("Log+Scale+PCA") + theme(plot.title = element_text(size = 10)) +
  scale_color_manual(values = c("A" = "#FF0000", # Bright red
                                "B" = "#0000FF", # Bright blue
                                "C" = "#CCCCCC", # Light grey
                                "D" = "#999999"))  # Darker grey
p_lpcascale <- p



### SCTRANSFORM (SEURAT)
colnames(Y) <- 1:J
rownames(Y) <- 1:I
Sco <- CreateSeuratObject(counts=Y)
Sco$group <- true_cluster
Sco <- SCTransform(Sco)
Sco <- RunPCA(Sco)
sct <- Sco@reductions$pca@cell.embeddings
df <- data.frame(x=sct[,1],y=sct[,2],color=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=pt.size)
p <- p + theme_bw()+xlab("")+ylab("")+guides(color="none") +
  ggtitle("SCT+PCA") +
  scale_color_manual(values = c("A" = "#FF0000", # Bright red
                                "B" = "#0000FF", # Bright blue
                                "C" = "#CCCCCC", # Light grey
                                "D" = "#999999"))  # Darker grey
p_sct <- p


### ANALYTIC PEARSON RESIDUALS
apr <- sctransform::vst(Y, method="offset")
pca.apr <- prcomp(t(apr$y))
df <- data.frame(x=pca.apr$x[,1],y=pca.apr$x[,2],color=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=pt.size)
p <- p + theme_bw()+xlab("")+ylab("")+guides(color="none")+
  ggtitle("APR+PCA") +
  scale_color_manual(values = c("A" = "#FF0000", # Bright red
                                "B" = "#0000FF", # Bright blue
                                "C" = "#CCCCCC", # Light grey
                                "D" = "#999999"))  # Darker grey
p_apr <- p

### LOG +  PCA
colnames(Y) <- 1:J
rownames(Y) <- 1:I
YL <- log(sweep(Y,MARGIN=2,STATS=L^{-1}*colSums(Y),FUN="/") + 1)
my.pca <- prcomp(t(YL))
lpca <- my.pca$x
df <- data.frame(x=lpca[,1],y=lpca[,2],color=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=pt.size)
p <- p + theme_bw()+xlab("")+ylab("")+guides(color="none") +
  ggtitle("Log+PCA") +
  scale_color_manual(values = c("A" = "#FF0000", # Bright red
                               "B" = "#0000FF", # Bright blue
                               "C" = "#CCCCCC", # Light grey
                               "D" = "#999999"))  # Darker grey
p_lpca <- p


library(ggpubr)
p.single.full <- ggarrange(ggarrange(pg1, pg2, pg3, nrow=1),
          ggarrange(p_lpca, p_lpcascale,
                    p_sct, p_apr, nrow=1), nrow=2,
        heights=c(1,1))
ggsave(p.single.full,
       filename="../plots/single_marker_pca.png",
       width=11.9, height=6.32, units="in")

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

