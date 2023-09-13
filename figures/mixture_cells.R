
##Where to save data
data.dir <- 0 #PUT YOUR DIRECTORY HERE

set.seed(1)
library(Seurat)
library(DuoClustering2018)
library(scGBM)

sce <- sce_full_Zhengmix8eq()

Y <- sce@assays@data@listData$counts
ct1 <- "b.cells"
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
ixs <- which(rowSums(Y.new) > 10)
Y.new <- Y.new[ixs,]

I <- nrow(Y.new)
J <- ncol(Y.new)
colnames(Y.new) <- 1:J
Sco <- CreateSeuratObject(counts=Y.new)
Sco <- NormalizeData(Sco)
Sco <- FindVariableFeatures(Sco)
Sco <- ScaleData(Sco)
Sco <- RunPCA(Sco)
Sco$my.type <- my.type
lpca <- Sco@reductions$pca@cell.embeddings

Sco <- RunUMAP(Sco, dims=1:20)
lpca.umap <- Sco@reductions$umap@cell.embeddings



Sco <- SCTransform(Sco)
Sco <- RunPCA(Sco,assay="SCT")
sct <- Sco@reductions$pca@cell.embeddings
df <- data.frame(x=sct[,1],y=sct[,2],color=my.type)
Sco <- RunUMAP(Sco, dims=1:20)
sct.umap <- Sco@reductions$umap@cell.embeddings


#Y.new <- Sco@assays$RNA@counts[Sco@assays$RNA@var.features,]
Y.new <- as.matrix(Y.new)
out <- gbm.sc(Y.new,M=20,max.iter=100, infer.beta = TRUE)

cell.mean <- colMeans(Y.new)

saveRDS(out, file=paste0(data.dir,"scGBM_gradient.RData"))

library(ggplot2)

###r2 analysis
gbm.r21 <- cor(my.type,out$scores[,1])^2
gbm.r22 <- cor(my.type,out$scores[,2])^2
sct.r21 <- cor(my.type,sct.umap[,1])^2
sct.r22 <- cor(my.type,sct.umap[,2])^2
l.r21 <- cor(my.type,lpca.umap[,1])^2
l.r22 <- cor(my.type,lpca.umap[,2])^2

lpca.r21 <- cor(my.type,lpca[,1])^2
lpca.r22 <- cor(my.type,lpca[,2])^2
sctpca.r21 <- cor(my.type,sct[,1])^2
sctpca.r22 <- cor(my.type,sct[,2])^2


##Plotting Figure

out <- readRDS(file = paste0(data.dir,"scGBM_gradient.RData"))

pt.size <- 0.5

df <- data.frame(x=sct.umap[,1],y=sct.umap[,2],color=my.type)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=pt.size)
p  <- p + scale_color_gradientn(colors=rainbow(2))
p <- p + theme_bw()+ guides(color="none")
p <- p + theme_bw()+xlab(paste0("UMAP1 (r2 = ", round(sct.r21,digits=2), ")"))
p <- p + ylab(paste0("UMAP2 (r2 = ", round(sct.r22,digits=2), ")"))
p <- p +  ggtitle("SCT+PCA+UMAP")
p_sct <- p


df <- data.frame(x=sct[,1],y=sct[,2],color=my.type)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=pt.size)
p  <- p + scale_color_gradientn(colors=rainbow(2))
p <- p + theme_bw()+ guides(color="none")
p <- p + theme_bw()+xlab(paste0("PCA1 (r2 = ", round(sctpca.r21,digits=2), ")"))
p <- p + ylab(paste0("PCA2 (r2 = ", round(sctpca.r22,digits=2), ")"))
p <- p +  ggtitle("SCT+PCA")
p_sctPCA <- p

df <- data.frame(x=my.type,y=sct[,1],color=my.type)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(alpha=0.5, size=pt.size)
p  <- p + scale_color_gradientn(colors=rainbow(2))
p <- p + theme_bw() + guides(color="none")+xlab("Mixture proportion")+ylab("")
p_sct.f1 <- p

df <- data.frame(x=log(cell.mean),y=sct[,1])
p <- ggplot(data=df,aes(x=x,y=y))+geom_point(alpha=0.5,color="grey",size=pt.size)
p <- p + theme_bw() + xlab("Log mean count")+ylab("")+ggtitle("SCT+PCA")
p <- p + stat_smooth(method="lm", formula=y~x, geom="smooth",
                     linetype="dashed", se=F)
p_sct.cm <- p

df <- data.frame(x=lpca.umap[,1],y=lpca.umap[,2],color=my.type)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=pt.size)
p  <- p + scale_color_gradientn(colors=rainbow(2))
p <- p + theme_bw()+xlab(paste0("UMAP1 (r2 = ", round(l.r21,digits=2), ")"))
p <- p + ylab(paste0("UMAP2 (r2 = ", round(l.r22,digits=2), ")"))
p <- p +  ggtitle("log+scale+PCA+UMAP")
p <- p + labs(color="Naive T")
p_log2 <- p


df <- data.frame(x=lpca[,1],y=lpca[,2],color=my.type)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=pt.size)
p  <- p + scale_color_gradientn(colors=rainbow(2))
p <- p + theme_bw()+xlab(paste0("PCA1 (r2 = ", round(lpca.r21,digits=2), ")"))
p <- p + ylab(paste0("PCA2 (r2 = ", round(lpca.r22,digits=2), ")"))
p <- p +  ggtitle("log+scale+PCA")
p <- p + labs(color="Naive T")
p_log2PCA <- p


df <- data.frame(x=my.type,y=lpca[,1],color=my.type)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(alpha=0.5,size=pt.size)
p  <- p + scale_color_gradientn(colors=rainbow(2))
p <- p + theme_bw() + guides(color="none")+xlab("")+ylab("")
p_log2.f1 <- p

df <- data.frame(x=log(cell.mean),y=lpca[,1])
p <- ggplot(data=df,aes(x=x,y=y))+geom_point(alpha=0.5,color="grey",size=pt.size)
p <- p + stat_smooth(method="lm", formula=y~x, geom="smooth",
                     linetype="dashed", se=F)
p <- p + theme_bw() + xlab("")+ylab("")+ggtitle("log+scale+PCA")
p_lpca.cm <- p

df <- data.frame(x=out$scores[,1],y=out$scores[,2],color=my.type)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=pt.size)
p  <- p + scale_color_gradientn(colors=rainbow(2))
p <- p + theme_bw() + guides(color="none")
p <- p + theme_bw()+xlab(paste0("GBM1 (r2 = ", round(gbm.r21,digits=2), ")"))
p <- p + ylab(paste0("GBM2 (r2 = ", round(gbm.r22,digits=2), ")"))
p <- p +  ggtitle("scGBM")
p_gbm <- p

df <- data.frame(x=my.type,y=out$scores[,1],color=my.type)
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(alpha=0.5,size=pt.size)
p  <- p + scale_color_gradientn(colors=rainbow(2))
p <- p + theme_bw() + guides(color="none")
p <- p + xlab("")+ylab("Factor 1 score")
p_gbm.f1 <- p

df <- data.frame(x=log(cell.mean),y=out$scores[,1])
p <- ggplot(data=df,aes(x=x,y=y))+geom_point(alpha=0.5,color="grey",size=pt.size)
p <- p + stat_smooth(method="lm", formula=y~x, geom="smooth",
                     linetype="dashed", se=F)
p <- p + theme_bw() + xlab("")+ylab("Factor 1 score")+ggtitle("scGBM")
p_gbm.cm <- p



library(ggpubr)
#p <- ggarrange(p_gbm,p_sctPCA,p_log2PCA,
#               p_gbm,p_sct, p_log2
#               p_gbm.f1,p_sct.f1,p_log2.f1,
#               p_gbm.cm, p_sct.cm, p_lpca.cm,nrow=3,ncol=3,
#               labels=c("a","","","b","","","c","",""))





#Volcano plot
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rownames(out$loadings)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart=mart)

ix <- match(G_list$ensembl_gene_id, rownames(out$loadings))
rownames(out$loadings)[ix] <- G_list$hgnc_symbol
out <- get.se(out,EPS=0.001)
p.volc <- loadings.volcano(out, return.plot=TRUE)

library(ggpubr)
pmid <- ggarrange(ggarrange(p_gbm,p_sct,p_log2,nrow=1,ncol=3,
                            labels=c("a","",""), common.legend=TRUE,
                            legend="bottom"),
                  ggarrange(p_gbm.cm,p_sct.cm,p_lpca.cm,nrow=1,ncol=3,
                            labels=c("b","","")),nrow=2,ncol=1,
                  heights=c(6,4))

p <- ggarrange(p_gbm,p_sct,p_log2,
               p_gbm.cm, p_sct.cm, p_lpca.cm,
               nrow=2,ncol=3,
               labels=c("a","","","b","",""))

library(gridExtra)
p <- grid.arrange(pmid, p.volc, ncol=1, heights=c(1, 0.75))

ggsave(p,filename=paste0(data.dir,"gradient.pdf"), units="in", width=6.5,
       height=8)
