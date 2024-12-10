
library(tidyverse)
set.seed(1)
library(Seurat)
library(viridis)
pt.size <- 0.5

set.seed(1)
I <- 1000
J <- 1000
spike.size <- 50
#baseline.means <- 10^{seq(-2, 3,length.out=I)}
baseline.means <- rep(1, I)
Mu <- matrix(baseline.means,nrow=I,ncol=J)
Mu[1,] <- 1
Mu[1,1:spike.size] <- 10
Y <- matrix(rpois(n=I*J,lambda=as.vector(Mu)),nrow=I,ncol=J)
colnames(Y) <- 1:ncol(Y); rownames(Y) <- 1:nrow(Y)


true_cluster <- rep(1,J)
true_cluster[1:spike.size] <- "A"
true_cluster[(spike.size + 1):1000] <- "B"
true_cluster <- as.character(true_cluster)



### ANALYTIC PEARSON RESIDUALS
apr <- sctransform::vst(Y[rowSums(Y) >= 5,], method="offset")
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

library(scGBM)

out <- gbm.sc(Y[rowSums(Y) >= 5,],M=20,max.iter=250,tol=10^{-10})

df <- data.frame(x=out$scores[,1], y=out$scores[,2],color=true_cluster)

p <- ggplot(data=df,aes(x=x,y=y,color=color))+geom_point(size=pt.size)
p_sm <- p + theme_bw()+xlab("GBM1")+ylab("GBM1")+guides(color="none") +
  scale_color_manual(values = c("A" = "#FF0000", # Bright red
                                "B" = "#0000FF"))




fit <- lm(ifelse(true_cluster == "A", 1, 0) ~ ., data=pca.apr$x[,c(1:20)] |> as.data.frame())

fit <- lm(ifelse(true_cluster == "A", 1, 0) ~ ., data=out$scores[,c(1:2)] |> as.data.frame())


generate_spike <- function(I, J, baseline.mean, spike.size, spike.mean) {
  Mu <- matrix(baseline.mean,nrow=I,ncol=J)
  Mu[1,] <- baseline.mean
  Mu[1,1:spike.size] <- spike.mean
  Y <- matrix(rpois(n=I*J,lambda=as.vector(Mu)),nrow=I,ncol=J)
  colnames(Y) <- 1:ncol(Y); rownames(Y) <- 1:nrow(Y)
  return(Y)
}

test_scGBM <- function(Y) {
  M <- 20
  out <- gbm.sc(Y,M=M)
  my.cor <- rep(1:M)
  u.true <- rep(0, nrow(Y)); u.true[1] <- 1
  for(m in 1:M) {
    my.cor[m] <- cor(out$U[,m], u.true)^2
  }
  return(max(my.cor))
}

test_apr <- function(Y) {
  M <- 20
  apr <- sctransform::vst(Y, method="offset")
  pca.apr <- irlba::prcomp_irlba(t(apr$y), n=20)
  my.cor <- rep(1:M)
  u.true <- rep(0, nrow(Y)); u.true[1] <- 1
  for(m in 1:M) {
    my.cor[m] <- cor(pca.apr$rotation[,m], u.true)^2
  }

  return(max(my.cor))
}

test_sct <- function(Y) {
  M <- 20
  J <- ncol(Y); I <- nrow(Y)
  colnames(Y) <- 1:J
  rownames(Y) <- 1:I
  Sco <- CreateSeuratObject(counts=Y)
  Sco <- SCTransform(Sco)
  pca.sct <- irlba::prcomp_irlba(t(Sco@assays$SCT$scale.data), n=20)
  my.cor <- rep(1:M)
  u.true <- rep(0, nrow(Y)); u.true[1] <- 1
  for(m in 1:M) {
    my.cor[m] <- cor(pca.sct$rotation[,m], u.true)^2
  }

  return(max(my.cor))
}

test_logpca <- function(Y) {
  M <- 20
  J <- ncol(Y); I <- nrow(Y)
  colnames(Y) <- 1:J
  rownames(Y) <- 1:I
  L <- median(colSums(Y))
  YL <- log(sweep(Y,MARGIN=2,STATS=L^{-1}*colSums(Y),FUN="/") + 1)
  pca.log <- irlba::prcomp_irlba(t(YL), n=20)
  my.cor <- rep(1:M)
  u.true <- rep(0, nrow(Y)); u.true[1] <- 1
  for(m in 1:M) {
    my.cor[m] <- cor(pca.log$rotation[,m], u.true)^2
  }

  return(max(my.cor))
}

baseline.means <- c(1, 100) #Small and large 
spike.fc <- c(0.01, 0.1, 10)
#iter <- 1:10 #10 repitions
iter <- c(1)
spike.size <- c(10, 25, 50, 100)

params <- expand.grid(baseline.means,
                      spike.fc,
                      spike.size,
                      iter)

params$scgbm <- rep(0, nrow(params))
params$apr <- rep(0, nrow(params))
params$sct <- rep(0, nrow(params))
params$logpca <- rep(0, nrow(params))

for(i in 1:nrow(params)) {
  cat("ITERATION ", i, "\n")
  Y <- generate_spike(I = 1000,J=1000,
                      baseline.mean = params$Var1[i], 
                      spike.size = params$Var3[i], 
                      spike.mean = params$Var1[i] * params$Var2[i])
  params$scgbm[i] <- test_scGBM(Y)
  params$apr[i] <- test_apr(Y)
  params$sct[i] <- test_sct(Y)
  params$logpca[i] <- test_logpca(Y)
}