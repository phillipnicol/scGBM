


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

Sco <- NormalizeData(Sco)
Sco <- FindVariableFeatures(Sco)
Sco <- ScaleData(Sco)
Sco <- RunPCA(Sco)
lpca <- Sco@reductions$pca@cell.embeddings
df <- data.frame(x=lpca[,1],y=lpca[,2],color=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=points.size)
p <- p + theme_bw()+xlab("PC1")+ylab("PC2")+labs(color="Cell Type")
p <- p + ggtitle("log+PCA")
p_log2 <- p

U <- Sco@reductions$pca@feature.loadings[,1]
df <- data.frame(x=1:I,y=U,color=ifelse(1:I<3,"red","grey"))
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=points.size)
p <- p + scale_color_manual(values=c("grey","red"))
p <- p + theme_bw()+xlab("Gene")+ylab("Factor 1 weight")+guides(color="none")
p <- p + ggtitle("log+PCA")
p_log2U <- p

Sco <- SCTransform(Sco)
Sco <- RunPCA(Sco)
sct <- Sco@reductions$pca@cell.embeddings
df <- data.frame(x=sct[,1],y=sct[,2],color=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=points.size)
p <- p + theme_bw()+xlab("SCT1")+ylab("SCT2")+guides(color="none")
p <- p + ggtitle("scTransform")
p_sct <- p

U <- Sco@reductions$pca@feature.loadings[,1]
df <- data.frame(x=1:I,y=U,color=ifelse(1:I<3,"red","grey"))
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=points.size)
p <- p + scale_color_manual(values=c("grey","red"))
p <- p + theme_bw()+xlab("Gene")+ylab("Factor 1 weight")+guides(color="none")
p <- p + ggtitle("scTransform")
p_sctU <- p


out <- gbm.sc(Y,M=50,tol=10^{-10})

df <- data.frame(x=out$V[,1],y=out$V[,2],color=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=points.size)
p <- p + theme_bw() + guides(color="none")
p <- p+xlab("GBM1")+ylab("GBM2")+ggtitle("scGBM")
p_gbm <- p

U <- Sco@reductions$pca@feature.loadings[,1]
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


df <- data.frame(x=V[,1],y=V[,2],color=true_cluster)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=points.size)
p <- p + theme_bw() + guides(color="none")
p <- p+xlab("GBM1")+ylab("GBM2")+ggtitle("scGBM")
p_glmpca <- p

df <- data.frame(x=1:I,y=U[,1],color=ifelse(1:I<3,"red","grey"))
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=points.size)
p <- p + scale_color_manual(values=c("grey","red"))
p <- p + theme_bw()+xlab("Gene")+ylab("Loading")+guides(color="none")
p_glmpcaU <- p
#Uniform perturbation

set.seed(1)
I <- 1000
J <- 1000

mu <- rexp(n=I,rate=1/10)

Mu <- matrix(mu,nrow=I,ncol=J)
Mu[,501:1000] <- sample(c(1/2,2),size=I,replace=TRUE)*Mu[,501:1000]

Y <- matrix(rpois(n=I*J,lambda=Mu),nrow=I,ncol=J)

out <- gbm.sc(Y,M=50)
df <- data.frame(x=mu,y=out$U[,1])
p <- ggplot(data=df,aes(x=x,y=y))
p <- p + geom_point(size=points.size)
p <- p + xlab(expression(mu)) + ylab("Factor 1 weight")+theme_bw() +ggtitle("Simulated")
p_unif_perturb <- p



p <- ggarrange(ggarrange(p_gbm,p_sct,p_log2,nrow=1,ncol=3,labels=c("a","","")),
               ggarrange(p_gbmU,p_sctU,p_log2U,nrow=1,ncol=3,labels=c("b","","")),
               ggarrange(p_unif_perturb,nrow=1,ncol=1,labels=c("c","d")),
               nrow=3,ncol=1)

p


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
