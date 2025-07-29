
set.seed(42)
library(Seurat)
expr <- ReadMtx("../../data/Zheng_ERCC/matrix.mtx",
                cells="../../data/Zheng_ERCC/barcodes.tsv",
                features="../../data/Zheng_ERCC/genes.tsv")


nz <- apply(expr, 1, function(x) sum(x != 0))

expr <- expr[nz >= 5,]

row.multiplier <- rexp(n=nrow(expr), rate=0.1)

expr2 <- sweep(expr, 1, row.multiplier, `*`)
expr2 <- as.matrix(expr2)

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
umap.apr.scaled <- umap::umap(my.pca$x)
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
umap.apr.scaled <- umap::umap(my.pca$x)
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

ggarrange(p.apr, p.apr.scaled, p.umap.scaled,
          p.sct, p.sct.scaled, p.sct.umap.scaled,
          nrow=2,ncol=3)



## scGBM with prior on sigma
out <- gbm.sc(expr2,M=20,sigma=10)

p.gbm2 <- data.frame(x=out$scores[,1], y=out$scores[,2]) |>
  ggplot(aes(x=x,y=y)) +
  geom_point(size=0.5) +
  theme_bw() +
  xlab("GBM1") + ylab("GBM2") +
  ggtitle("ERCC Scaled (GBM)")

ggarrange(p.gbm1,p.gbm2,nrow=1)








library(tidyverse)
set.seed(1)
library(Seurat)
pt.size <- 0.5

set.seed(1)
I <- 1000
J <- 1000
Mu <- matrix(0.1,nrow=I,ncol=J)

true_cluster <- rep(0,J)
for(c in 1:10) {
  z <- rexp(n=1)
  up <- sample(c(1,10),size=2,replace=FALSE)
  Mu[c,] <- up[1]*z
  Mu[c,(100*(c-1) + 1):(100*c)] <- up[2]*z
  true_cluster[(100*(c-1) + 1):(100*c)] <- as.character(c)
}


#Mu[2,667:1000] <- 50
Y <- matrix(rpois(n=I*J,lambda=as.vector(Mu)),nrow=I,ncol=J)

true_cluster <- rep(1,J)
true_cluster[1:979] <- "C"
true_cluster[980:990] <- "A"
true_cluster[991:1000] <- "B"


df <- data.frame(x=out$scores[,1],y=out$scores[,2],color=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=pt.size)
p <- p + theme_bw()+xlab("SCT1")+ylab("SCT2")+guides(color="none")
p_sct <- p

rownames(Y) <- 1:1000; colnames(Y) <- 1:1000
apr <- sctransform::vst(Y, method="offset")
my.pca <- irlba::prcomp_irlba(apr$y,n=20)
umap.apr <- umap::umap(my.pca$rotation)
df <- data.frame(x=my.pca$rotation[,1],y=my.pca$rotation[,2],color=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=pt.size)
p <- p + theme_bw()+xlab("SCT1")+ylab("SCT2")+guides(color="none")
p_sct <- p

out <- gbm.sc(Y,M=10)
df <- data.frame(x=out$scores[,1],y=out$scores[,2],color=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=pt.size)
p <- p + theme_bw()

