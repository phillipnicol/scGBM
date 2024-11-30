
library(tidyverse)
set.seed(1)
library(Seurat)
library(viridis)
pt.size <- 0.5

set.seed(1)
I <- 1000
J <- 1000
spike.size <-50
#baseline.means <- rexp(n=I, rate=1)
baseline.means <- rep(1, I)
Mu <- matrix(baseline.means,nrow=I,ncol=J)
Mu[1,] <- 100
Mu[1,1:spike.size] <- 0
Y <- matrix(rpois(n=I*J,lambda=as.vector(Mu)),nrow=I,ncol=J)
colnames(Y) <- 1:ncol(Y); rownames(Y) <- 1:nrow(Y)


true_cluster <- rep(1,J)
true_cluster[1:spike.size] <- "A"
true_cluster[(spike.size + 1):1000] <- "B"
true_cluster <- as.character(true_cluster)



### ANALYTIC PEARSON RESIDUALS
apr <- sctransform::vst(Y, method="offset")
pca.apr <- prcomp(t(apr$y))
apr.umap <- umap::umap(pca.apr$x[,1:10])$layout
df <- data.frame(x=pca.apr$x[,1],y=pca.apr$x[,2],color=true_cluster) |> arrange(desc(color))
p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=pt.size)
p <- p + theme_bw()+xlab("")+ylab("")+guides(color="none")+
  ggtitle("APR+PCA") +
  scale_color_manual(values = c("A" = "#FF0000", # Bright red
                                "B" = "#0000FF") # Bright blue
  )
p_apr <- p


#devtools::install_github("phillipnicol/scGBM",ref="dev2")


out <- gbm.sc(Y,M=20)

df <- data.frame(x=out$scores[,3], y=out$scores[,4],color=true_cluster)

p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=pt.size)
p_sm <- p + theme_bw()+xlab("GBM1")+ylab("GBM1")+guides(color="none") +
  scale_color_manual(values = c("A" = "#FF0000", # Bright red
                                "B" = "#0000FF"))
