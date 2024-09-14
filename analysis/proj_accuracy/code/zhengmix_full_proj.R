
library(DuoClustering2018)
library(scGBM)

sce <- sce_full_Zhengmix8eq()

Y <- counts(sce)
Y <- as.matrix(Y)
set.seed(1)
J <- ncol(Y)
#perm <- sample(1:J,size=J,replace=FALSE)
#Rmse <- matrix(0,nrow=length(subset),ncol=iters)

library(fastglm)
library(scGBM)
library(rstiefel)
library(doParallel)

Y <- Y[rowSums(Y) >= 20, ]
out <- gbm.sc(Y,M=20)

out2 <- gbm.sc(Y,M=20, subset=3994,ncores=6)

my.cor <- rep(0,20)
for(m in 1:20) {
  cat("Factor ", m, " Cor: ",
      cor(out$scores[,m], out2$scores[,m])^2,
      " ", out$D[m], "\n")
  my.cor[m] <- abs(cor(out$scores[,m], out2$scores[,m]))
}

Proj.diff <- sqrt(sum(((out$V %*% t(out$V) - out2$scores %*% solve(t(out2$scores)%*% out2$scores) %*% t(out2$scores))^2)))

p <- data.frame(x=out$D[1:10],
                 y=my.cor[1:10]) |>
  ggplot(aes(x=x,y=y)) + geom_point() +
  theme_bw() + xlab(expression(sigma)) + ylab("Magnitude of correlation")

ggsave(p, filename="../plots/full_subsetsize_comparison.png")
