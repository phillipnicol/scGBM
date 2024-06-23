library(rstiefel)
library(glmpca)
library(scGBM)
library(fastglm)

I <- 10^3
J <- 10^4
M <- 10

reps <- 100

#Fisher, Avagrad, SGD, scGBM, scGBM-proj
nmethod <- 5

res <- array(0,dim=c(reps,nmethod,M)) #For factor cor
res2 <- matrix(0,nrow=reps, ncol=nmethod) #For X mse

find.best.cor <- function(V.true, factor.test) {
  val <- rep(0, ncol(V.true))
  for(m in 1:ncol(V.true)) {
    val[m] <- abs(cor(V.true[,m], factor.test))
  }
  return(max(val))
}

set.seed(1)
for(i in 1:reps) {
  sim <- simData(I,J,d=M)
  Pv.true <- sim$V %*% solve(t(sim$V) %*% sim$V) %*% t(sim$V)

  #gbm-sc
  out <- gbm.sc(sim$Y,M=10,order.by.deviance = FALSE)
  for(m in 1:M) {
    res[i,1,m] <- find.best.cor(sim$V,out$scores[,m])
  }
  Pv.pred <- out$V %*% t(out$V)
  res2[i,1] <- sqrt(mean((Pv.true-Pv.pred)^2))

  out <- gbm.sc(sim$Y,M=10,order.by.deviance = FALSE,
                subset=round(0.1*J), ncores=8)
  for(m in 1:M) {
    res[i,2,m] <- find.best.cor(sim$V,out$scores[,m])
  }
  Pv.pred <- out$scores %*% solve(t(out$scores) %*% out$scores) %*% t(out$scores)
  res2[i,2] <- sqrt(mean((Pv.true-Pv.pred)^2))

  fit <- glmpca(Y=sim$Y,L=10)
  fit$res$factors <- as.matrix(fit$res$factors)
  fit$res$loadings <- as.matrix(fit$res$loadings)
  for(m in 1:M) {
    res[i,3,m] <- find.best.cor(sim$V,fit$res$factors[,m])
  }
  Pv.pred <- fit$res$factors %*% solve(t(fit$res$factors) %*% fit$res$factors) %*% t(fit$res$factors)
  X.pred <- fit$res$loadings %*% t(fit$res$factors)
  res2[i,3] <- sqrt(mean((Pv.true-Pv.pred)^2))

  print("SGD")
  print(typeof(sim$Y))
  rownames(sim$Y) <- 1:I; colnames(sim$Y) <- 1:J
  fit <- glmpca(Y=sim$Y,Y.oos=sim$Y,L=10,minibatch="stochastic",ctl=list(batch_size=0.1*J))
  fit$res$factors <- as.matrix(fit$res$factors)
  fit$res$loadings <- as.matrix(fit$res$loadings)
  for(m in 1:M) {
    res[i,4,m] <- find.best.cor(sim$V,fit$res$factors[,m])
  }
  Pv.pred <- fit$res$factors %*% solve(t(fit$res$factors) %*% fit$res$factors) %*% t(fit$res$factors)
  X.pred <- fit$res$loadings %*% t(fit$res$factors)
  res2[i,4] <-sqrt(mean((Pv.true-Pv.pred)^2))

  fit <- glmpca(Y=sim$Y, L=10,optimizer = "fisher")
  fit$res$factors <- as.matrix(fit$res$factors)
  fit$res$loadings <- as.matrix(fit$res$loadings)
  for(m in 1:M) {
    res[i,5,m] <- find.best.cor(sim$V,fit$res$factors[,m])
  }
  X.pred <- fit$res$loadings %*% t(fit$res$factors)
  Pv.pred <- fit$res$factors %*% solve(t(fit$res$factors) %*% fit$res$factors) %*% t(fit$res$factors)
  res2[i,5] <- sqrt(mean((Pv.true-Pv.pred)^2))
}

saveRDS(res, "../data/factor_correlation.RDS")
saveRDS(res2, "../data/mse.RDS")


### Plotting
res <- readRDS("../data/factor_correlation.RDS")
res2 <- readRDS("../data/mse.RDS")

df <- reshape2::melt(res)

library(tidyverse)
method.names = c("scGBM-full",
                 "scGBM-proj",
                 "GLM-PCA (AvaGrad)",
                 "GLM-PCA (SGD)",
                 "GLM-PCA (Fisher)")

df <- df |> mutate(Method = method.names[Var2]) |>
  group_by(Method, Var3) |>
  summarise(mean=mean(value^2),
            ymin=mean(value^2) - sd(value^2),
            ymax=mean(value^2) + sd(value^2))

p <- ggplot(df, aes(x=Var3, y=mean, color=Method,ymin=ymin,ymax=ymax)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  xlab("Latent factor") +
  ylab("r^2 with ground truth")


df <- reshape2::melt(res2)
df$value <- sqrt(I*J)*df$value
df <- df |> mutate(Method = method.names[Var2])
p <- ggplot(df,aes(x=Method, y=value, fill=Method)) +
  geom_boxplot(alpha=0.5) +
  geom_jitter(shape=16,position=position_jitter(0.2)) +
  ylab(expression(paste("||", Pi[hat(V)], " - ", Pi[V], "||"))) +
  theme_bw() + xlab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  guides(fill="none")

ggsave(p, filename="../plots/sim_accuracy.png")
