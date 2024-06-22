

library(Seurat)
library(tidyverse)
set.seed(1)
expr <- ReadMtx("../../data/Zheng_ERCC/matrix.mtx",
                cells="../../data/Zheng_ERCC/barcodes.tsv",
                features="../../data/Zheng_ERCC/genes.tsv")


nz <- apply(expr, 1, function(x) sum(x != 0))

expr <- expr[nz >= 5,]

Y <- as.matrix(expr)

gene.means <- rep(0, nrow(Y))
gene.vars <- rep(0, nrow(Y))

for(i in 1:nrow(Y)) {
  gene.means[i] <- mean(mean(colSums(Y))*Y[i,]/colSums(Y))
  gene.vars[i] <- var(mean(colSums(Y))*Y[i,]/colSums(Y))
}

plain <- function(x,...) {
  format(x, ..., scientific = FALSE, drop0trailing=TRUE)
}

df <- data.frame(x=gene.means,y=gene.vars) |>
  ggplot(aes(x=x,y=y)) +
  geom_point(size=1) +
  theme_bw() +
  xlab("Gene mean") +
  ylab("Gene variance") +
  scale_x_log10(labels=plain) +
  scale_y_log10(labels=plain)
  #ggtitle("Mean variance relationship: ERCC controls")

ggsave(df, filename="../plots/mean_variance_relationship.png")

## Default ScTransform
apr <- sctransform::vst(expr)
my.pca <- irlba::prcomp_irlba(apr$y)
umap.sct <- umap::umap(my.pca$rotation)
p.sct <- data.frame(x=my.pca$x[,1], y=my.pca$x[,2]) |>
  ggplot(aes(x=x,y=y)) +
  geom_point(size=0.5) +
  xlab("PCA1") + ylab("PCA2") +
  ggtitle("ERCC (SCT + PCA)") +
  theme_bw()


## APR
apr <- sctransform::vst(expr, method="offset")
my.pca <- irlba::prcomp_irlba(apr$y,n=10)
umap.apr <- umap::umap(my.pca$rotation)
p.apr <- data.frame(x=my.pca$x[,1], y=my.pca$x[,2]) |>
  ggplot(aes(x=x,y=y)) +
  geom_point(size=0.5) +
  xlab("PCA1") + ylab("PCA2") +
  ggtitle("ERCC (APR + PCA)") +
  theme_bw()

## Log + scale + PCA
size.factor <- colSums(expr)
CPM <- 10^6*sweep(expr, 2, size.factor, "/")
log2CPM <- log2(CPM+1)
my.pca <- prcomp(log2CPM)
#umap.log.scale <- umap::umap(my.pca$rotation[,1:10])

## Log + PCA
size.factor <- colSums(expr)
normalized <- sweep(expr, 2, size.factor, "/")
log2.counts <- log(normalized + 1,base=2)
my.pca <- irlba::prcomp_irlba(log2.counts)
umap.log2 <- umap::umap(my.pca$rotation)

p <- data.frame(x=log(colSums(expr)), y=my.pca$rotation[,1],
                color=log(colSums(expr))) |>
  ggplot(aes(x=x,y=y,color=color)) + geom_point() +
  scale_color_gradient(low="blue", high="red")

## scGBM with prior on sigma
out <- gbm.sc(expr |> as.matrix(),M=20,sigma=10)




### Make it even worse by multiplying each row
#set.seed(1)
row.multiplier <- rexp(n=nrow(expr), rate=0.1)
#row.multiplier <- rep(100, nrow(expr))

expr2 <- sweep(expr, 1, row.multiplier, `*`)
expr2 <- as.matrix(expr2)


## Default ScTransform
apr <- sctransform::vst(expr2)
my.pca <- irlba::prcomp_irlba(t(apr$y))
umap.sct <- umap::umap(my.pca$x)
p.sct.scaled <- data.frame(x=my.pca$x[,1], y=my.pca$x[,2]) |>
  ggplot(aes(x=x,y=y)) +
  geom_point(size=0.5) +
  xlab("PCA1") + ylab("PCA2") +
  ggtitle("ERCC scaled (SCT + PCA)") +
  theme_bw()


## APR
apr <- sctransform::vst(expr2, method="offset")
my.pca <- irlba::prcomp_irlba(t(apr$y))
umap.apr <- umap::umap(my.pca$x)
p.apr.scaled <- data.frame(x=my.pca$x[,1], y=my.pca$x[,2]) |>
  ggplot(aes(x=x,y=y)) +
  geom_point(size=0.5) +
  xlab("PCA1") + ylab("PCA2") +
  ggtitle("ERCC scaled (APR + PCA)") +
  theme_bw()

## Log + PCA
size.factor <- colSums(expr2)
normalized <- sweep(expr2, 2, size.factor, "/")
log2.counts <- log(normalized + 1,base=2)
my.pca <- irlba::prcomp_irlba(log2.counts)
umap.log2 <- umap::umap(my.pca$rotation)

## Log + scale + PCA
size.factor <- colSums(expr2)
CPM <- 10^6*sweep(expr2, 2, size.factor, "/")
log2CPM <- log2(CPM+1)
my.pca <- prcomp(log2CPM)

## scGBM with prior on sigma
out.2 <- gbm.sc(expr2,M=10,sigma=10)



library(ggpubr)

p.ercc.full <- ggarrange(p.sct, p.sct.scaled,
               p.apr, p.apr.scaled,
               nrow=2, ncol=2)


ggsave(p.ercc.full,filename="../plots/ercc_full.png")



## Make scGBM plot

p.gbm.unscaled <- data.frame(x=out$scores[,1], y=out$scores[,2]) |>
  ggplot(aes(x=x,y=y)) +
  geom_point(size=0.5) +
  xlab("GBM1") + ylab("GBM2") +
  ggtitle("ERCC (GBM)") +
  theme_bw()


p.gbm.scaled <- data.frame(x=out.2$scores[,1], y=out.2$scores[,2]) |>
  ggplot(aes(x=x,y=y)) +
  geom_point(size=0.5) +
  xlab("GBM1") + ylab("GBM2") +
  ggtitle("ERCC scaled (GBM)") +
  theme_bw()

p.ercc.gbm <- ggarrange(p.gbm.unscaled,p.gbm.scaled,
                         nrow=1, ncol=2)

ggsave(p.ercc.gbm,filename="../plots/ercc_gbm.png",
       width=8.72, height=3.47, units="in")
