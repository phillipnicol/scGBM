

library(Seurat)
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


## Default ScTransform
apr <- sctransform::vst(expr)
my.pca <- irlba::prcomp_irlba(apr$y)
umap.sct <- umap::umap(my.pca$rotation)


## APR
apr <- sctransform::vst(expr, method="offset")
my.pca <- irlba::prcomp_irlba(apr$y,n=10)
umap.apr <- umap::umap(my.pca$rotation)

## Log + scale + PCA
size.factor <- colSums(expr)
CPM <- 10^6*sweep(expr2, 2, size.factor, "/")
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
set.seed(1)
row.multiplier <- rexp(n=nrow(expr), rate=0.01)

expr2 <- sweep(expr, 1, row.multiplier, `*`)
expr2 <- as.matrix(expr2)


## Default ScTransform
apr <- sctransform::vst(expr2)
my.pca <- irlba::prcomp_irlba(t(apr$y))
umap.sct <- umap::umap(my.pca$x)


## APR
apr <- sctransform::vst(expr2, method="offset")
my.pca <- irlba::prcomp_irlba(t(apr$y))
umap.apr <- umap::umap(my.pca$rotation)

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
out <- gbm.sc(expr2,M=10,sigma=1)













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

out <- gbm.sc(Y,M=2)
fit <- fastglmpca::fit_glmpca_pois(Y,K=20)

