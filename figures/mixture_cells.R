

### Mixture cell example

set.seed(1)
library(Seurat)
library(DuoClustering2018)

sce <- sce_full_Zhengmix8eq()

Y <- sce@assays$data@listData$counts

ct1 <- "memory.t"
ct2 <- "naive.t"
Y <- Y[,sce$phenoid %in% c(ct1,ct2)]
cell_type <- sce$phenoid[sce$phenoid %in% c(ct1,ct2)]

ixs <- apply(Y,1,function(x) sum(x!=0))
ixs <- which(ixs > 10)

Y <- Y[ixs,]

cell_type <- sce$phenoid[sce$phenoid %in% c(ct1,ct2)]


ct1.ix <- which(cell_type == ct1)
ct2.ix <- which(cell_type==ct2)

J <- 5000
Y.new <- matrix(0,nrow=nrow(Y),ncol=J)
rownames(Y.new) <- rownames(Y)
c <- 1
my.type <- rep(0,J)
for(j in 1:J) {
  counts.new <- rep(0,nrow(Y))
  ix <- sample(1:nrow(Y),size=1)
  cell <- sample(ct1.ix,size=1)
  counts.new <- Y[,cell]
  cell <- sample(ct2.ix,size=1)
  gene <- sample(1:nrow(Y),size=ix,replace=FALSE)
  counts.new[gene] <- Y[gene,cell]
  Y.new[,j] <- counts.new
  my.type[j] <- ix/nrow(Y)
}


I <- nrow(Y.new)
J <- ncol(Y.new)
colnames(Y.new) <- 1:J
Sco <- CreateSeuratObject(counts=Y.new)
Sco <- NormalizeData(Sco)
Sco <- FindVariableFeatures(Sco,nfeatures = 6088)
Sco <- ScaleData(Sco)
Sco <- RunPCA(Sco)
Sco$my.type <- my.type
lpca <- Sco@reductions$pca@cell.embeddings


Sco <- SCTransform(Sco,
                   variable.features.n = 6088)
Sco <- RunPCA(Sco,assay="SCT")
sct <- Sco@reductions$pca@cell.embeddings
df <- data.frame(x=sct[,1],y=sct[,2],color=my.type)


Y.new <- Sco@assays$RNA@counts[Sco@assays$RNA@var.features,]
ixs <- which(rowSums(Y.new) > 10)
Y.new <- Y.new[ixs,]
Y.new <- as.matrix(Y.new)
out <- gbm.sc(Y.new,M=20,max.iter =50)

cell.mean <- colMeans(Y.new)





##Plotting Figure

df <- data.frame(x=sct[,1],y=sct[,2],color=my.type)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point()
p  <- p + scale_color_gradientn(colors=rainbow(2))
p <- p + theme_bw()+ guides(color="none")+xlab("SCT1")+ylab("SCT2")+ggtitle("SCTransform")
p_sct <- p

df <- data.frame(x=my.type,y=sct[,1],color=my.type)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(alpha=0.5)
p  <- p + scale_color_gradientn(colors=rainbow(2))
p <- p + theme_bw() + guides(color="none")+xlab("Mixture proportion")+ylab("")
p_sct.f1 <- p

df <- data.frame(x=log(cell.mean),y=sct[,1])
p <- ggplot(data=df,aes(x=x,y=y))+geom_point(alpha=0.5,color="grey")
p <- p + theme_bw() + xlab("Log mean count")+ylab("")
p_sct.cm <- p

df <- data.frame(x=lpca[,1],y=lpca[,2],color=my.type)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point()
p  <- p + scale_color_gradientn(colors=rainbow(2))
p <- p + theme_bw()+xlab("PC1")+ylab("PC2")+ggtitle("Log2+PCA")
p <- p + labs(color="Naive T")
p_log2 <- p

df <- data.frame(x=my.type,y=lpca[,1],color=my.type)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(alpha=0.5)
p  <- p + scale_color_gradientn(colors=rainbow(2))
p <- p + theme_bw() + guides(color="none")+xlab("")+ylab("")
p_log2.f1 <- p

df <- data.frame(x=log(cell.mean),y=lpca[,1])
p <- ggplot(data=df,aes(x=x,y=y))+geom_point(alpha=0.5,color="grey")
p <- p + theme_bw() + xlab("")+ylab("")
p_lpca.cm <- p

df <- data.frame(x=out$V[,1],y=out$V[,2],color=my.type)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point()
p  <- p + scale_color_gradientn(colors=rainbow(2))
p <- p + theme_bw() + guides(color="none")
p <- p+xlab("GBM1")+ylab("GBM2")+ggtitle("scGBM")
p_gbm <- p

df <- data.frame(x=my.type,y=out$V[,1],color=my.type)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(alpha=0.5)
p  <- p + scale_color_gradientn(colors=rainbow(2))
p <- p + theme_bw() + guides(color="none")
p <- p + xlab("")+ylab("Component 1 score")
p_gbm.f1 <- p

df <- data.frame(x=log(cell.mean),y=out$V[,1])
p <- ggplot(data=df,aes(x=x,y=y))+geom_point(alpha=0.5,color="grey")
p <- p + theme_bw() + xlab("")+ylab("Component 1 score")
p_gbm.cm <- p


library(ggpubr)
p <- ggarrange(p_gbm,p_sct,p_log2,
               p_gbm.f1,p_sct.f1,p_log2.f1,
               p_gbm.cm, p_sct.cm, p_lpca.cm,nrow=3,ncol=3,
               labels=c("a","","","b","","","c","",""))
p
