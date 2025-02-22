
library(tidyverse)
set.seed(1)
library(Seurat)
library(viridis)
pt.size <- 0.5
library(scGBM)

generate_spike <- function(I, J, baseline.mean, spike.size, spike.mean, second.spike.mean) {
  Mu <- matrix(1,nrow=I,ncol=J)
  Mu[1,] <- baseline.mean
  Mu[1,1:spike.size] <- spike.mean
  Mu[2,] <- 1
  Mu[2,667:1000] <- second.spike.mean
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
    my.cor[m] <- out$U[1,m]/max(out$U[-1,m])
    #my.cor[m] <- cor(out$U[,m], u.true)^2
  }
  my.cor <- abs(my.cor)

  return(max(my.cor))
}

test_apr <- function(Y) {
  M <- 20
  apr <- sctransform::vst(Y, method="offset")
  pca.apr <- irlba::prcomp_irlba(t(apr$y), n=20)
  my.cor <- rep(1:M)
  u.true <- rep(0, nrow(Y)); u.true[1] <- 1
  for(m in 1:M) {
    my.cor[m] <- pca.apr$rotation[1,m]/max(pca.apr$rotation[-1,m])
    #my.cor[m] <- cor(pca.apr$rotation[,m], u.true)^2
  }

  my.cor <- abs(my.cor)

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
    my.cor[m] <- pca.sct$rotation[1,m]/max(pca.sct$rotation[-1,m])
    #my.cor[m] <- cor(pca.sct$rotation[,m], u.true)^2
  }

  my.cor <- abs(my.cor)

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
    my.cor[m] <- pca.log$rotation[1,m]/max(pca.log$rotation[-1,m])
    #my.cor[m] <- cor(pca.log$rotation[,m], u.true)^2
  }

  my.cor <- abs(my.cor)

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
    my.cor[m] <- pca.log$rotation[1,m]/max(pca.log$rotation[-1,m])
    #my.cor[m] <- cor(pca.log$rotation[,m], u.true)^2
  }

  my.cor <- abs(my.cor)

  return(max(my.cor))
}

test_logpca_S1 <- function(Y) {
  M <- 20
  J <- ncol(Y); I <- nrow(Y)
  colnames(Y) <- 1:J
  rownames(Y) <- 1:I
  #L <- median(colSums(Y))
  YL <- log(sweep(Y,MARGIN=2,STATS=L^{-1},FUN="/") + 1)
  pca.log <- irlba::prcomp_irlba(t(YL), n=20)
  my.cor <- rep(1:M)
  u.true <- rep(0, nrow(Y)); u.true[1] <- 1
  for(m in 1:M) {
    my.cor[m] <- pca.log$rotation[1,m]/max(pca.log$rotation[-1,m])
    #my.cor[m] <- cor(pca.log$rotation[,m], u.true)^2
  }

  my.cor <- abs(my.cor)

  return(max(my.cor))
}

baseline.means <- 100 #Small and large
spike.mean <- c(10, 20, 50)
#iter <- 1:10 #10 repitions
iter <- c(1:10)
spike.size <- c(1,10, 25, 50, 100)
second.spike.mean <- c(10, 50, 100)

params <- expand.grid(baseline.means,
                      spike.mean,
                      spike.size,
                      iter,
                      second.spike.mean)

params$scgbm <- rep(0, nrow(params))
params$apr <- rep(0, nrow(params))
params$sct <- rep(0, nrow(params))
params$logpca <- rep(0, nrow(params))

for(i in 1:nrow(params)) {
  cat("ITERATION ", i, "\n")
  Y <- generate_spike(I = 1000,J=1000,
                      baseline.mean=params$Var1[i],
                      spike.size=params$Var3[i],
                      spike.mean=params$Var2[i],
                      second.spike.mean=params$Var5[i])
  params$scgbm[i] <- test_scGBM(Y)
  params$apr[i] <- test_apr(Y)
  params$sct[i] <- test_sct(Y)
  params$logpca[i] <- test_logpca(Y)
}

saveRDS(params,file="../data/spike_in_params.RDS")




### Plotting

params <- readRDS("../data/spike_in_params.RDS")

params <- params[,-1] #Remove baseline means

colnames(params)[1:4] <- c("spike.mean",
                      "spike.size",
                      "replicate",
                      "second.spike.mean")

library(ggplot2)
library(reshape2)

colnames(params)[5:8] <- c("scGBM", "APR+PCA", "SCT+PCA","Log+PCA")

df <- reshape2::melt(params,measure.vars=c("scGBM", "APR+PCA", "SCT+PCA","Log+PCA"))

#colnames(df) <- c("Marker mean", "# of cell type A", "Replicate", "Second.spike.mean",
#                  "Method", "value")

df <- df |> group_by(spike.mean,spike.size,variable,second.spike.mean) |>
  summarize(mean=mean(value)) |>
  ggplot(aes(x=spike.size, y=mean,color=variable)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  xlab("# of cell type A") + ylab("Separation") +
  labs(color = "Method") +
  xlim(c(0,50)) + ylim(c(0,1.1)) +
  geom_abline(slope=0, intercept=1, color="grey", linetype="dashed") +
  facet_grid(second.spike.mean ~ spike.mean)

ggsave(df, filename="../plots/spike_in_sim.png",
       width=8.67, height=5.57, units="in")



