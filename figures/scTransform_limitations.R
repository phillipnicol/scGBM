


library(tidyverse)
set.seed(1)
library(Seurat)
pt.size <- 0.5

set.seed(1)
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

Sco <- SCTransform(Sco)
Sco <- RunPCA(Sco)
sct <- Sco@reductions$pca@cell.embeddings
df <- data.frame(x=sct[,1],y=sct[,2],color=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=pt.size)
p <- p + theme_bw()+xlab("SCT1")+ylab("SCT2")+guides(color="none")
p_sct <- p


res <- Sco@assays$SCT@scale.data
rownames(res)[1:2] <- c("Gene 1", "Gene 2")
rownames(res)[3:1000] <- "Null genes"
res <- t(res)
df <- reshape2::melt(res)
df$Var2 <- factor(df$Var2,levels=c("Null genes", "Gene 2", "Gene 1"))
p <- ggplot(data=df,aes(x=Var2,y=value,fill=Var2))
p <- p + geom_boxplot()+guides(fill="none")
p <- p + scale_fill_manual(values=c("lavender","lightsalmon","forestgreen"))
p <- p + xlab("")+ylab("Pearson residual")
p <- p + theme_bw() + coord_flip()
p_sctgp <- p


U <- Sco@reductions$pca@feature.loadings[,1]
df <- data.frame(x=1:I,y=U,color=ifelse(1:I<3,"red","grey"))
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=pt.size)
p <- p + scale_color_manual(values=c("grey","red"))
p <- p + theme_bw()+xlab("Gene")+ylab("PC1 Weight")+guides(color="none")
p_sctU <- p

df <- data.frame(y=Y[1,],x=true_cluster,fill=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,fill=fill))+geom_boxplot()
p <- p + theme_bw() + guides(fill="none")
p <- p + xlab("") + ylab("Counts") + ggtitle("Gene 1")
p_g1 <- p

df <- data.frame(y=Y[2,],x=true_cluster,fill=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,fill=fill))+geom_boxplot()
p <- p + theme_bw()+ guides(fill="none")
p <- p + xlab("Cell Type") + ylab("Counts") + ggtitle("Gene 2")
p_g2 <- p

df <- data.frame(y=Y[3,],x=true_cluster,fill=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,fill=fill))+geom_boxplot()
p <- p + theme_bw()+guides(fill="none")
p <- p + xlab("") + ylab("Counts") + ggtitle("Genes 3-1000 (Random noise)")
p_g3 <- p



#Uniform perturbation

set.seed(1)
I <- 1000
J <- 1000

mu <- rexp(n=I,rate=1/10)

Mu <- matrix(mu,nrow=I,ncol=J)
Mu[,501:1000] <- sample(c(1/2,2),size=I,replace=TRUE)*Mu[,501:1000]

Y <- matrix(rpois(n=I*J,lambda=Mu),nrow=I,ncol=J)


rownames(Y) <- 1:I; colnames(Y) <- 1:J

library(Seurat)

Sco <- CreateSeuratObject(counts=Y)
Sco <- SCTransform(Sco)
Sco <- RunPCA(Sco)
loadings1 <- Sco@reductions$pca@feature.loadings[rownames(Y),1]


df <- data.frame(mu=mu,sct=loadings1)

library(ggplot2)

p <- ggplot(data=df,aes(x=mu,y=loadings1))
p <- p + geom_point(size=pt.size)
p <- p + xlab(expression(mu)) + ylab("PC1 Weight") + theme_bw() +ggtitle("Simulation (ii)")
p_unif_perturb <- p


library(ggpubr)
p <- ggarrange(ggarrange(p_g1,p_g2,p_g3,nrow=1,ncol=3,labels=c("a","","")),
               ggarrange(p_sct,p_sctU,p_sctgp,nrow=1,ncol=3,labels=c("b","c","d")),
               ggarrange(p_unif_perturb,nrow=1,ncol=1,labels=c("e","f")),
               nrow=3)



### For the example on Monocytes
### Download genes.tsv and barcodes.tsv from
### the 10X website

library(Seurat)
genes <- read.csv("genes.tsv")
barcodes <- read.csv("barcodes.tsv")

expr <- ReadMtx("matrix.mtx",cells="barcodes.tsv",
                features="genes.tsv")

Sco <-Seurat::CreateSeuratObject(counts=expr)

Sco <- SCTransform(Sco)
Sco <- RunPCA(Sco)

Y <- Sco@assays$RNA@counts[Sco@assays$SCT@var.features,]
#CPM <-  sweep(Y, 2, 10^6/colSums(Y), "*")
gene.means <- rowMeans(Y)
df <- data.frame(x=log(gene.means),y=Sco@reductions$pca@feature.loadings[,1])
p <- ggplot(data=df,aes(x=x,y=y))
p <- p + geom_point(size=pt.size)
p <- p + xlab("Log Mean Count") + ylab("PC1 Weight") + ggtitle("Monocytes")+theme_bw()
p_sct_monocytes <- p
