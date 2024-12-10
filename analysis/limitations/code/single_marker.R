
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
Mu[1,1:10] <- 20
Mu[1,11:20] <- 10
Mu[2,] <- 1
Mu[2,667:1000] <- 2000
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
                               "C" = "#FFD700", # Light grey
                               "D" = "#999999"))+  # Darker grey
  scale_y_sqrt()

df <- data.frame(y=Y[2,],x=true_cluster,fill=true_cluster)
pg2 <- ggplot(data=df,aes(x=x,y=y,fill=fill))+geom_boxplot()+
  theme_bw() + guides(fill="none") +
  geom_jitter(alpha=0.5, size=0.5,width=0.3) +
  xlab("") + ylab("Counts") + ggtitle("Gene 2") +
  scale_fill_manual(values = c("A" = "#FF0000", # Bright red
                                "B" = "#0000FF", # Bright blue
                                "C" = "#FFD700", # Light grey
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
                               "C" = "#FFD700", # Light grey
                               "D" = "#999999"))+  # Darker grey
  scale_y_sqrt()


### LOG + SCALE + PCA (SEURAT)
colnames(Y) <- 1:J
rownames(Y) <- 1:I
Sco <- CreateSeuratObject(counts=Y)
Sco <- NormalizeData(Sco)
Sco <- FindVariableFeatures(Sco)
Sco <- ScaleData(Sco,scale=TRUE)
Sco$group <- true_cluster
Sco <- RunPCA(Sco)
lpca <- Sco@reductions$pca@cell.embeddings
lpca.scale.umap <- umap::umap(lpca[,1:10])$layout
df <- data.frame(x=lpca[,1],y=lpca[,2],color=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=pt.size)
p <- p + theme_bw()+xlab("")+ylab("")+guides(color="none") +
  ggtitle("Log+Scale+PCA") + theme(plot.title = element_text(size = 10)) +
  scale_color_manual(values = c("A" = "#FF0000", # Bright red
                                "B" = "#0000FF", # Bright blue
                                "C" = "#FFD700", # Light grey
                                "D" = "#999999"))  # Darker grey
p_lpcascale <- p

pe_lpcascale <- data.frame(x=2:50, y = Sco@reductions$pca@stdev[-1]) |>
  ggplot(aes(x=x,y=y)) + geom_point() +
  theme_bw() +
  geom_hline(yintercept = 1, color="red", linetype="dashed") +
  xlab("PC") +
  ylab(expression(sqrt(lambda))) +
  ggtitle("Log+Scale+PCA")

p_lpscale_umap <- data.frame(x=lpca.scale.umap[,1], y=lpca.scale.umap[,2],color=true_cluster) |>
  ggplot(aes(x=x,y=y,color=color))+geom_point(size=pt.size) +
  theme_bw()+xlab("")+ylab("")+guides(color="none") +
  ggtitle("Log+Scale+PCA+UMAP") + theme(plot.title = element_text(size = 10)) +
  scale_color_manual(values = c("A" = "#FF0000", # Bright red
                                "B" = "#0000FF", # Bright blue
                                "C" = "#FFD700", # Light grey
                                "D" = "#999999"))  # Darker grey


### SCTRANSFORM (SEURAT)
colnames(Y) <- 1:J
rownames(Y) <- 1:I
Sco <- CreateSeuratObject(counts=Y)
Sco$group <- true_cluster
Sco <- SCTransform(Sco)
Sco <- RunPCA(Sco)
sct <- Sco@reductions$pca@cell.embeddings
sct.umap <- umap::umap(sct[,1:10])$layout
df <- data.frame(x=sct[,1],y=sct[,2],color=true_cluster) |> arrange(desc(color))
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=pt.size)
p <- p + theme_bw()+xlab("")+ylab("")+guides(color="none") +
  ggtitle("SCT+PCA") +
  scale_color_manual(values = c("A" = "#FF0000", # Bright red
                                "B" = "#0000FF", # Bright blue
                                "C" = "#FFD700", # Light grey
                                "D" = "#999999"))  # Darker grey
p_sct <- p

pe_sct <- data.frame(x=2:50, y = Sco@reductions$pca@stdev[-1]) |>
  ggplot(aes(x=x,y=y)) + geom_point() +
  theme_bw() +
  geom_hline(yintercept = sd(Sco@assays$SCT@scale.data[1,]), color="red", linetype="dashed") +
  xlab("PC") +
  ylab(expression(sqrt(lambda))) +
  ggtitle("SCT+PCA")


p_sct_umap <- data.frame(x=sct.umap[,1], y=sct.umap[,2],color=true_cluster) |>
  ggplot(aes(x=x,y=y,color=color))+geom_point(size=pt.size) +
  theme_bw()+xlab("")+ylab("")+guides(color="none") +
  ggtitle("SCT+PCA+UMAP") + theme(plot.title = element_text(size = 10)) +
  scale_color_manual(values = c("A" = "#FF0000", # Bright red
                                "B" = "#0000FF", # Bright blue
                                "C" = "#FFD700", # Light grey
                                "D" = "#999999"))  # Darker grey



### ANALYTIC PEARSON RESIDUALS
apr <- sctransform::vst(Y, method="offset")
pca.apr <- prcomp(t(apr$y))
apr.umap <- umap::umap(pca.apr$x[,1:10])$layout
df <- data.frame(x=pca.apr$x[,1],y=pca.apr$x[,2],color=true_cluster) |> arrange(desc(color))
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=pt.size)
p <- p + theme_bw()+xlab("")+ylab("")+guides(color="none")+
  ggtitle("APR+PCA") +
  scale_color_manual(values = c("A" = "#FF0000", # Bright red
                                "B" = "#0000FF", # Bright blue
                                "C" = "#FFD700", # Light grey
                                "D" = "#999999"))  # Darker grey
p_apr <- p

pe_apr <- data.frame(x=2:50, y = pca.apr$sdev[2:50]) |>
  ggplot(aes(x=x,y=y)) + geom_point() +
  theme_bw() +
  geom_hline(yintercept = sd(apr$y[1,]), color="red", linetype="dashed") +
  xlab("PC") +
  ylab(expression(sqrt(lambda))) +
  ggtitle("APR+PCA")

p_apr_umap <- data.frame(x=apr.umap[,1], y=apr.umap[,2],color=true_cluster) |>
  ggplot(aes(x=x,y=y,color=color))+geom_point(size=pt.size) +
  theme_bw()+xlab("")+ylab("")+guides(color="none") +
  ggtitle("APR+PCA+UMAP") + theme(plot.title = element_text(size = 10)) +
  scale_color_manual(values = c("A" = "#FF0000", # Bright red
                                "B" = "#0000FF", # Bright blue
                                "C" = "#FFD700", # Light grey
                                "D" = "#999999"))  # Darker grey



### LOG +  PCA
colnames(Y) <- 1:J
rownames(Y) <- 1:I
L <- median(colSums(Y))
YL <- log(sweep(Y,MARGIN=2,STATS=L^{-1}*colSums(Y),FUN="/") + 1)
my.pca <- prcomp(t(YL))
lpca <- my.pca$x
lpca.umap <- umap::umap(lpca[,1:10])$layout
df <- data.frame(x=lpca[,1],y=lpca[,2],color=true_cluster) |> arrange(desc(color))
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=pt.size)
p <- p + theme_bw()+xlab("")+ylab("")+guides(color="none") +
  ggtitle("Log+PCA") +
  scale_color_manual(values = c("A" = "#FF0000", # Bright red
                               "B" = "#0000FF", # Bright blue
                               "C" = "#FFD700", # Light grey
                               "D" = "#999999"))  # Darker grey
p_lpca <- p

pe_lpca <- data.frame(x=2:50, y = my.pca$sdev[2:50]) |>
  ggplot(aes(x=x,y=y)) + geom_point() +
  theme_bw() +
  geom_hline(yintercept = sd(YL[1,]), color="red", linetype="dashed") +
  xlab("PC") +
  ylab(expression(sqrt(lambda))) +
  ggtitle("Log+PCA")

p_lpca_umap <- data.frame(x=lpca.umap[,1], y=lpca.umap[,2],color=true_cluster) |>
  ggplot(aes(x=x,y=y,color=color))+geom_point(size=pt.size) +
  theme_bw()+xlab("")+ylab("")+guides(color="none") +
  ggtitle("LOG+PCA+UMAP") + theme(plot.title = element_text(size = 10)) +
  scale_color_manual(values = c("A" = "#FF0000", # Bright red
                                "B" = "#0000FF", # Bright blue
                                "C" = "#FFD700", # Light grey
                                "D" = "#999999"))  # Darker grey



library(ggpubr)
p.single.full <- ggarrange(ggarrange(pg1, pg2, pg3, nrow=1),
          ggarrange(p_lpca, p_lpcascale,
                    p_sct, p_apr, nrow=1), nrow=2,
        heights=c(1,1), labels=c("a","b"))
ggsave(p.single.full,
       filename="../plots/single_marker_pca.png",
       width=11.9, height=6.32, units="in")

p.pe <- ggarrange(pe_lpca, pe_lpcascale,
                  pe_sct, pe_apr, ncol=2, nrow=2)

ggsave(p.pe,
       filename="../plots/single_marker_eigenvalues.png",
       width=8.5,height=8.17, units="in")


p.umap <- ggarrange(p_lpca_umap, p_lpscale_umap,
                    p_sct_umap, p_apr_umap, nrow=2,ncol=2)


ggsave(p.umap,
       filename="../plots/single_marker_umap.png",
       width=8.5,height=8.17, units="in")




set.seed(42)
library(Seurat)
expr <- ReadMtx("../../data/Zheng_ERCC/matrix.mtx",
                cells="../../data/Zheng_ERCC/barcodes.tsv",
                features="../../data/Zheng_ERCC/genes.tsv")


nz <- apply(expr, 1, function(x) sum(x != 0))

expr <- expr[nz >= 5,]

row.multiplier <- rexp(n=nrow(expr), rate=0.1)

#expr2 <- sweep(expr, 1, row.multiplier, `*`)
#expr2 <- as.matrix(expr2)

expr2 <- expr
for(i in 1:nrow(expr2)) {
  expr2[i,] <- expr2[i,]/mean(expr2[i,])
}

## APR + PCA + (UMAP)

apr <- sctransform::vst(expr, method="offset")
my.pca <- irlba::prcomp_irlba(t(apr$y),n=10)
umap.apr <- umap::umap(my.pca$x)
p.apr <- data.frame(x=my.pca$x[,1], y=my.pca$x[,2]) |>
  ggplot(aes(x=x,y=y)) +
  geom_point(size=0.5) +
  xlab("PCA1") + ylab("PCA2") +
  ggtitle("ERCC (APR + PCA)") +
  theme_bw() +
  xlab("") + ylab("")

apr <- sctransform::vst(expr2, method="offset")
my.pca.scaled <- irlba::prcomp_irlba(t(apr$y),n=10)
umap.apr.scaled <- umap::umap(my.pca.scaled$x)
p.apr.scaled <- data.frame(x=my.pca.scaled$x[,1], y=my.pca.scaled$x[,2]) |>
  ggplot(aes(x=x,y=y)) +
  geom_point(size=0.5) +
  xlab("PCA1") + ylab("PCA2") +
  ggtitle("ERCC Scaled (APR + PCA)") +
  theme_bw() +
  xlab("") + ylab("")

p.umap.scaled <- data.frame(x=umap.apr.scaled$layout[,1],
                            y=umap.apr.scaled$layout[,2]) |>
  ggplot(aes(x=x,y=y)) +
  geom_point(size=0.5) +
  xlab("PCA1") + ylab("PCA2") +
  ggtitle("ERCC Scaled (APR + PCA + UMAP)") +
  theme_bw() +
  xlab("") + ylab("")

## SCT + PCA + (UMAP)

apr <- sctransform::vst(expr)
my.pca <- irlba::prcomp_irlba(t(apr$y),n=10)
umap.apr <- umap::umap(my.pca$x)
p.sct <- data.frame(x=my.pca$x[,1], y=my.pca$x[,2]) |>
  ggplot(aes(x=x,y=y)) +
  geom_point(size=0.5) +
  xlab("PCA1") + ylab("PCA2") +
  ggtitle("ERCC (SCT + PCA)") +
  theme_bw() +
  xlab("") + ylab("")

apr <- sctransform::vst(expr2)
my.pca.scaled <- irlba::prcomp_irlba(t(apr$y),n=10)
umap.apr.scaled <- umap::umap(my.pca.scaled$x)
p.sct.scaled <- data.frame(x=my.pca.scaled$x[,1], y=my.pca.scaled$x[,2]) |>
  ggplot(aes(x=x,y=y)) +
  geom_point(size=0.5) +
  xlab("PCA1") + ylab("PCA2") +
  ggtitle("ERCC Scaled (SCT + PCA)") +
  theme_bw() +
  xlab("") + ylab("")

p.sct.umap.scaled <- data.frame(x=umap.apr.scaled$layout[,1],
                                y=umap.apr.scaled$layout[,2]) |>
  ggplot(aes(x=x,y=y)) +
  geom_point(size=0.5) +
  xlab("PCA1") + ylab("PCA2") +
  ggtitle("ERCC Scaled (SCT + PCA + UMAP)") +
  theme_bw() +
  xlab("") + ylab("")


library(ggpubr)

p.ercc.scaled <- ggarrange(p.apr, p.apr.scaled, p.umap.scaled,
          p.sct, p.sct.scaled, p.sct.umap.scaled,
          nrow=2,ncol=3, labels="c")

library(ggpubr)
p <- ggarrange(p.single.full, p.ercc.scaled, nrow=2,
               heights=c(1,1))

ggsave(p,filename="../plots/limitation_plot.png",
       width=12.8, height=11.7)


###scGBM
set.seed(42)
out <- gbm.sc(Y,M=20,sigma=10, infer.beta=TRUE)


df <- data.frame(x=out$scores[,1], y=out$scores[,2],color=true_cluster)

p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=pt.size)
p_sm <- p + theme_bw()+xlab("GBM1")+ylab("GBM1")+guides(color="none") +
  scale_color_manual(values = c("A" = "#FF0000", # Bright red
                              "B" = "#0000FF", # Bright blue
                              "C" = "#FFD700", # Light grey
                              "D" = "#999999"))


## scGBM with prior on sigma
out <- gbm.sc(expr |> as.matrix(),M=20,sigma=10)

p.gbm1 <- data.frame(x=out$scores[,1], y=out$scores[,2]) |>
  ggplot(aes(x=x,y=y)) +
  geom_point(size=0.5) +
  theme_bw() +
  xlab("GBM1") + ylab("GBM2") +
  ggtitle("ERCC (GBM)")

out <- gbm.sc(expr2 |> as.matrix(),M=20,sigma=10)

p.gbm2 <- data.frame(x=out$scores[,1], y=out$scores[,2]) |>
  ggplot(aes(x=x,y=y)) +
  geom_point(size=0.5) +
  theme_bw() +
  xlab("GBM1") + ylab("GBM2") +
  ggtitle("ERCC Scaled (GBM)")

p_ercc <- ggarrange(p.gbm1,p.gbm2,nrow=1)

res <- readRDS("../../cluster_accuracy/data/zhengmix8eq.RDS")
res2 <- readRDS("../../cluster_accuracy/data/zhengmix8uneq.RDS")

names(res) <- c("scGBM-full", "scGBM-proj",
                "log+scale+PCA",
                "SCT+PCA",
                "APR+PCA",
                "GLM-PCA")

df <- data.frame(equal = res,
                 unequal = res2, Method=names(res) |> fct_inorder())

p <- reshape2::melt(df,id.vars="Method") |>
  ggplot(aes(x = Method, y=value, fill=variable)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(labels=c("Balanced", "Unbalanced"),
                    values=c("firebrick", "forestgreen")) +
  ylab("ARI") +
  labs(fill = "Size distribution") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(p, filename="../plots/gbm_limitations_excerpt_c.png",
       width=6.4, height=3.57,units="in")

p <- ggarrange(p_sm, p_ercc, p, nrow=3, labels=c("a","b","c"))

ggsave(p, filename="../plots/gbm_limitations.png",
       width=2541, height=3508, units="px")

ggsave(p.ercc.scaled, filename="../plots/limitation_plot_panel_c.png",
       width=12.8, height=5.5)

ggsave(p_ercc, filename = "../plots/gbm_limitations_excerpt_b.png",
       width=6.37, height=4.79, units="in")




