

load("../data/zhengmix_Time.RData")
load("../data/zhengmix_LL.RData")


library(tidyverse)
## Plotting function for loglike vs runtime comparison


#cutoff <- min(which(abs(rel.diff) < 10^{-4}))

df.1 <- data.frame(time=Time[[1]],
                   Method="scGBM-full",
                   ll=LL[[1]])


df.2 <- data.frame(time=rowMeans(Time[[2]]),
                   Method="scGBM-proj",
                   ll=rowMeans(LL[[2]]))


df.3 <- data.frame(time=Time[[3]],
                   Method="GLM-PCA (AvaGrad)",
                   ll=LL[[3]])

df.4 <- data.frame(time=Time[[4]],
                   Method="GLM-PCA (Fisher)",
                   ll=LL[[4]])

df.5 <- data.frame(time=Time[[5]],
                   Method="GLM-PCA (SGD)",
                   ll=LL[[5]])



df <- rbind(df.1,df.2,df.3,df.4,df.5) %>% as.data.frame


p <- ggplot(data=df,aes(x=time/3600,y=ll,color=Method))
p <- p + geom_point() + geom_line()
#p <- p + geom_segment(x=log10(Time[[2]][1]/3600),xend=100,y=LL[[2]][1],yend=LL[[2]][1],
#color=hue_pal()(5)[5],linetype="dashed")
p <- p + scale_x_log10(
  breaks = scales::trans_breaks("log10", function(x) 10^x),
  labels = scales::trans_format("log10", scales::math_format(10^.x)),
  limits=c(10,10^{3})/3600
)
p <- p + xlab("") + ylab("")
p <- p + theme_bw()
p <- p + annotation_logticks(sides = 'b')
p <- p + ggtitle("10X Immune (J=3,994; I=6,049)")
p <- p + guides(color="none")
p  <- p
p_tenximmune <- p



library(DuoClustering2018)
sce <- sce_full_Zhengmix8eq()


##Compare embedding plots
point.size <- 0.5


#scGBM
df <- readRDS("../data/gbm_zhengmix_embedding.RDS")[,c(1,2)]
colnames(df) <- c("x","y")
df <- df |> as.data.frame() |> mutate(color=sce$phenoid)
p_scgbm <- ggplot(data=df,aes(x=x,y=y,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("scGBM")

df <- readRDS("../data/gbm_zhengmix_umap.RDS")
df <- df |> as.data.frame() |> mutate(color=sce$phenoid)
p_scgbm_umap <- ggplot(data=df,aes(x=V1,y=V2,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("scGBM+UMAP")

#scGBM proj

df <- readRDS("../data/../data/gbm_proj_zhengmix_embedding.RDS") |> as.data.frame() |>
  mutate(color=sce$phenoid)
p_scgbm_proj <- ggplot(data=df,aes(x=V1,y=V2,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("scGBM-proj")

df <- readRDS("../data/../data/gbm_proj_zhengmix_umap.RDS") |> as.data.frame() |>
  mutate(color=sce$phenoid)
p_scgbm_proj_umap <- ggplot(data=df,aes(x=V1,y=V2,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("scGBM-proj+UMAP")

#Fisher
df <- readRDS("../data/glmpca_fisher_zhengmix_embedding.RDS")[,c(1,2)] |> as.data.frame() |>
  mutate(color=sce$phenoid)
p_glmpca_fisher <- ggplot(data=df,aes(x=dim1,y=dim2,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("GLM-PCA (Fisher)")

df <- readRDS("../data/glmpca_fisher_zhengmix_umap.RDS") |> as.data.frame() |>
  mutate(color=sce$phenoid)
p_glmpca_fisher_umap <- ggplot(data=df,aes(x=V1,y=V2,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("GLM-PCA (Fisher) + UMAP")

#Avagrad
df <- readRDS("../data/glmpca_avagrad_zhengmix_embedding.RDS")[,c(1,2)] |> as.data.frame() |>
  mutate(color=sce$phenoid)
p_glmpca_avagrad <- ggplot(data=df,aes(x=dim1,y=dim2,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("GLM-PCA (AvaGrad)")

df <- readRDS("../data/glmpca_avagrad_zhengmix_umap.RDS") |> as.data.frame() |>
  mutate(color=sce$phenoid)
p_glmpca_avagrad_umap <- ggplot(data=df,aes(x=V1,y=V2,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("GLM-PCA (Avagrad) + UMAP")


#SGD

df <- readRDS("../data/glmpca_sgd_zhengmix_embedding.RDS")[,c(1,2)] |> as.data.frame() |>
  mutate(color=sce$phenoid)
p_glmpca_sgd <- ggplot(data=df,aes(x=dim1,y=dim2,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("GLM-PCA (SGD)")


df <- readRDS("../data/glmpca_sgd_zhengmix_umap.RDS") |> as.data.frame() |>
  mutate(color=sce$phenoid)
p_glmpca_sgd_umap <- ggplot(data=df,aes(x=V1,y=V2,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("GLM-PCA (SGD) + UMAP")

df <- readRDS("../data/apr_zhengmix_embedding.RDS")[,c(1,2)] |> as.data.frame() |>
  mutate(color=sce$phenoid)
p_apr <- ggplot(data=df,aes(x=PC1,y=PC2,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("APR+PCA")

df <- readRDS("../data/apr_zhengmix_umap.RDS")[,c(1,2)] |> as.data.frame() |>
  mutate(color=sce$phenoid)
p_apr_umap <- ggplot(data=df,aes(x=V1,y=V2,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("APR+PCA+UMAP")

df <- readRDS("../data/sct_zhengmix_embedding.RDS")[,c(1,2)] |> as.data.frame() |>
  mutate(color=sce$phenoid)
p_sct <- ggplot(data=df,aes(x=PC_1,y=PC_2,color=color)) +
  geom_point(size=point.size) + theme_bw() + labs(color="Cell Type") +
  ggtitle("SCT+PCA") + xlab("") + ylab("")

df <- readRDS("../data/sct_zhengmix_umap.RDS")[,c(1,2)] |> as.data.frame() |>
  mutate(color=sce$phenoid)
p_sct_umap <- ggplot(data=df,aes(x=V1,y=V2,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("SCT+PCA+UMAP")


df <- readRDS("../data/logpca_zhemgmix_embedding.RDS")[,c(1,2)] |> as.data.frame() |>
  mutate(color=sce$phenoid)
p_lpca <- ggplot(data=df,aes(x=PC1,y=PC2,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("Log+PCA")


df <- readRDS("../data/logpca_zhengmix_umap.RDS")[,c(1,2)] |> as.data.frame() |>
  mutate(color=sce$phenoid)
p_lpca_umap <- ggplot(data=df,aes(x=V1,y=V2,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("Log+PCA+UMAP")


df <- readRDS("../data/seurat_zhengmix_embedding.RDS")[,c(1,2)] |> as.data.frame() |>
  mutate(color=sce$phenoid)
p_log_scale <- ggplot(data=df,aes(x=PC_1,y=PC_2,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("Log+Scale+PCA")

df <- readRDS("../data/seurat_zhengmix_umap.RDS")[,c(1,2)] |> as.data.frame() |>
  mutate(color=sce$phenoid)
p_log_scale_umap <- ggplot(data=df,aes(x=V1,y=V2,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("Log+Scale+PCA+UMAP")


library(ggpubr)

p.embed <- ggarrange(p_scgbm, p_scgbm_proj,
                     p_glmpca_avagrad, p_glmpca_fisher,
                     p_glmpca_sgd, p_apr,
                     p_sct, p_log_scale,
                     p_lpca, ncol=3, nrow=3,
                     common.legend=TRUE,
                     legend="bottom")

ggsave(p.embed, filename="../plots/10x_embeddings.png",
       width=10.4, height=8.36, units="in")

p.umap <- ggarrange(p_scgbm_umap, p_scgbm_proj_umap,
                    p_glmpca_avagrad_umap, p_glmpca_fisher_umap,
                    p_glmpca_sgd_umap, p_apr_umap,
                    p_sct_umap, p_log_scale_umap,
                    p_lpca_umap, ncol=3, nrow=3,
                    common.legend=TRUE,
                    legend="bottom")

ggsave(p.umap, filename="../plots/10x_umap.png",
       width=10.4, height=8.36, units="in")

#ggsave(p, filename="../plots/tenx_immune_embeddings.png",
#       width=1.5*10.6, height=1.5*8.33, units="in")



load("../data/blish_time.RData")
load("../data/blish_LL.RData")


library(tidyverse)
## Plotting function for loglike vs runtime comparison

ll <- LL[[1]]

df.1 <- data.frame(time=Time[[1]],
                   Method="scGBM-full",
                   ll=LL[[1]])

df.2 <- data.frame(time=Time[[2]],
                   Method="scGBM-proj",
                   ll=LL[[2]])



df.3 <- data.frame(time=Time[[3]],
                   Method="GLM-PCA (AvaGrad)",
                   ll=LL[[3]])


df.4 <- data.frame(time=Time[[4]],
                   Method="GLM-PCA (Fisher)",
                   ll=LL[[4]])

df.5 <- data.frame(time=Time[[5]],
                   Method="GLM-PCA (SGD)",
                   ll=LL[[5]])



df <- rbind(df.1,df.2,df.3,df.4, df.5) %>% as.data.frame


p <- ggplot(data=df,aes(x=time/3600,y=ll,color=Method))
p <- p + geom_point() + geom_line()
p <- p + scale_x_log10(
  breaks = scales::trans_breaks("log10", function(x) 10^x),
  labels = scales::trans_format("log10", scales::math_format(10^.x)),
  limits=c(10^2,10^{5})/3600
)
p <- p + xlab("") + ylab("") +
  ylim(-5*10^7, NA)
p <- p + annotation_logticks(sides = 'b')
p <- p + theme_bw()
p <- p + ggtitle("COVID-19 (J=44,721; I=17,393)")
pblish <- p


## Near zero result


df.1 <- data.frame(time=Time[[1]],
                   Method="scGBM-full",
                   ll=LL[[1]])

df.2 <- data.frame(time = near.zero.results$time.blish,
                   Method="scGBM-full (GLM-PCA init)",
                   ll=near.zero.results$ll.blish)


df <- rbind(df.1,df.2) %>% as.data.frame


p <- ggplot(data=df,aes(x=time/3600,y=ll,color=Method))
p <- p + geom_point() + geom_line()
#p <- p + geom_segment(x=log10(Time[[2]][1]/3600),xend=100,y=LL[[2]][1],yend=LL[[2]][1],
#color=hue_pal()(5)[5],linetype="dashed")
p <- p + scale_x_log10(
  breaks = scales::trans_breaks("log10", function(x) 10^x),
  labels = scales::trans_format("log10", scales::math_format(10^.x)),
  limits=c(10,10^{3})/3600
)
p <- p + xlab("") + ylab("")
p <- p + theme_bw()
p <- p + annotation_logticks(sides = 'b')
p <- p + ggtitle("10X Immune (J=3,994; I=6,049)")
#p <- p + guides(color="none")
p  <- p
p_tenximmune.initcompare <- p





##Compare embedding plots
blish_meta <- readRDS("../../data/blish/blish_meta.RDS")
point.size <- 0.33
df <- readRDS("../data/apr_blish_embedding.RDS")[,c(1,2)] |> as.data.frame() |>
  mutate(color=blish_meta$cell.type.coarse)
p_apr <- ggplot(data=df,aes(x=PC1,y=PC2,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("APR+PCA")

df <- readRDS("../data/apr_blish_umap.RDS")[,c(1,2)] |> as.data.frame() |>
  mutate(color=blish_meta$cell.type.coarse)
p_apr_umap <- ggplot(data=df,aes(x=V1,y=V2,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("APR+PCA+UMAP")

df <- readRDS("../data/sct_blish_embedding.RDS")[,c(1,2)] |> as.data.frame() |>
  mutate(color=blish_meta$cell.type.coarse)
p_sct <- ggplot(data=df,aes(x=PC_1,y=PC_2,color=color)) +
  geom_point(size=point.size) + theme_bw() + labs(color="Cell Type") +
  ggtitle("SCT+PCA")

df <- readRDS("../data/sct_blish_umap.RDS")[,c(1,2)] |> as.data.frame() |>
  mutate(color=blish_meta$cell.type.coarse)
p_sct_umap <- ggplot(data=df,aes(x=V1,y=V2,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("SCT+PCA+UMAP")


df <- readRDS("../data/logpca_blish_embedding.RDS")[,c(1,2)] |> as.data.frame() |>
  mutate(color=blish_meta$cell.type.coarse)
p_lpca <- ggplot(data=df,aes(x=PC1,y=PC2,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("Log+PCA")


df <- readRDS("../data/logpca_blish_umap.RDS")[,c(1,2)] |> as.data.frame() |>
  mutate(color=blish_meta$cell.type.coarse)
p_lpca_umap <- ggplot(data=df,aes(x=V1,y=V2,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("Log+PCA+UMAP")



#scGBM
df <- readRDS("../data/gbm_blish_embedding.RDS")[,c(1,2)]
colnames(df) <- c("x","y")
df <- df |> as.data.frame() |> mutate(color=blish_meta$cell.type.coarse)
p_scgbm <- ggplot(data=df,aes(x=x,y=y,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("scGBM")


df <- readRDS("../data/gbm_blish_umap.RDS")
df <- df |> as.data.frame() |> mutate(color=blish_meta$cell.type.coarse)
p_scgbm_umap <- ggplot(data=df,aes(x=V1,y=V2,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("scGBM+UMAP")

## scGBM-proj

df <- readRDS("../data/gbm_proj_blish_embedding.RDS") |> as.data.frame() |>
  mutate(color=blish_meta$cell.type.coarse)
p_scgbm_proj <- ggplot(data=df,aes(x=V1,y=V2,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("scGBM-proj")

df <- readRDS("../data/gbm_proj_blish_umap.RDS")
df <- df |> as.data.frame() |> mutate(color=blish_meta$cell.type.coarse)
p_scgbm_proj_umap <- ggplot(data=df,aes(x=V1,y=V2,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("scGBM-proj+UMAP")


#Avagrad

df <- readRDS("../data/glmpca_avagrad_blish_embedding.RDS")[,c(1,2)] |> as.data.frame() |>
  mutate(color=blish_meta$cell.type.coarse)
p_glmpca_avagrad <- ggplot(data=df,aes(x=dim1,y=dim2,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("GLM-PCA (AvaGrad)")

df <- readRDS("../data/glmpca_avagrad_blish_umap.RDS") |> as.data.frame() |>
  mutate(color=blish_meta$cell.type.coarse)
p_glmpca_avagrad_umap <- ggplot(data=df,aes(x=V1,y=V2,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("GLM-PCA (AvaGrad) + UMAP")

#SGD

df <- readRDS("../data/glmpca_sgd_blish_embedding.RDS")[,c(1,2)] |> as.data.frame() |>
  mutate(color=blish_meta$cell.type.coarse)
p_glmpca_sgd <- ggplot(data=df,aes(x=dim1,y=dim2,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("GLM-PCA (AvaGrad)")

df <- readRDS("../data/glmpca_sgd_blish_umap.RDS") |> as.data.frame() |>
  mutate(color=blish_meta$cell.type.coarse)
p_glmpca_sgd_umap <- ggplot(data=df,aes(x=V1,y=V2,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("GLM-PCA (SGD) + UMAP")

#Fisher

df <- readRDS("../data/glmpca_fisher_blish_embedding.RDS")[,c(1,2)] |> as.data.frame() |>
  mutate(color=blish_meta$cell.type.coarse)
p_glmpca_fisher <- ggplot(data=df,aes(x=dim1,y=dim2,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("GLM-PCA (Fisher)")

df <- readRDS("../data/glmpca_fisher_blish_umap.RDS") |> as.data.frame() |>
  mutate(color=blish_meta$cell.type.coarse)
p_glmpca_fisher_umap <- ggplot(data=df,aes(x=V1,y=V2,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("GLM-PCA (Fisher) + UMAP")

#SEURAT
df <- readRDS("../data/seurat_blish_embedding.RDS")[,c(1,2)] |> as.data.frame() |>
  mutate(color=blish_meta$cell.type.coarse)
p_log_scale <- ggplot(data=df,aes(x=PC_1,y=PC_2,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("Log+Scale+PCA")

df <- readRDS("../data/seurat_blish_umap.RDS")[,c(1,2)] |> as.data.frame() |>
  mutate(color=blish_meta$cell.type.coarse)
p_log_scale_umap <- ggplot(data=df,aes(x=V1,y=V2,color=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  ggtitle("Log+Scale+PCA+UMAP")

library(ggpubr)

p.embed <- ggarrange(p_scgbm, p_scgbm_proj,
                     p_glmpca_avagrad, p_glmpca_fisher,
                     p_glmpca_sgd, p_apr,
                     p_sct, p_log_scale,
                     p_lpca, ncol=3, nrow=3,
                     common.legend=TRUE,
                     legend="bottom")

ggsave(p.embed, filename="../plots/blish_embeddings.png",
       width=10.4, height=8.36, units="in")

p.umap <- ggarrange(p_scgbm_umap, p_scgbm_proj_umap,
                    p_glmpca_avagrad_umap, p_glmpca_fisher_umap,
                    p_glmpca_sgd_umap, p_apr_umap,
                    p_sct_umap, p_log_scale_umap,
                    p_lpca_umap, ncol=3, nrow=3,
                    common.legend=TRUE,
                    legend="bottom")

ggsave(p.umap, filename="../plots/blish_umap.png",
       width=10.4, height=8.36, units="in")

p <- ggarrange(p_scgbm, p_scgbm_proj, p_scgbm_umap, p_glmpca_avagrad,
               p_glmpca_fisher, p_glmpca_sgd,
               p_sct_umap, p_apr_umap, p_lpca_umap,
               nrow=3, ncol=3, common.legend = TRUE)

#ggsave(p, filename="../plots/blish_embeddings.png", width=12.22, height=9.5, units="in")


df <- readRDS("../data/apr_blish_embedding.RDS")[,c(1,2)] |> as.data.frame() |>
  mutate(color=ifelse(blish_meta$cell.type.coarse == "RBC", "RBC", "Other"))
p_apr_2 <- ggplot(data=df,aes(x=PC1,y=PC2,color=color,alpha=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  scale_color_manual(values=c("grey", "red"), labels=c("Other", "RBC")) +
  scale_alpha_manual(values=c(0.3, 1), labels=c("Other", "RBC")) +
  guides(alpha="none") +
  ggtitle("APR+PCA")

df <- readRDS("../data/apr_blish_umap.RDS")[,c(1,2)] |> as.data.frame() |>
  mutate(color=ifelse(blish_meta$cell.type.coarse == "RBC", "RBC", "Other"))
p_apr_umap_2 <- ggplot(data=df,aes(x=V1,y=V2,color=color,alpha=color)) +
  geom_point(size=point.size) + theme_bw() +
  xlab("") + ylab("") + labs(color="Cell Type") +
  scale_color_manual(values=c("grey", "red"), labels=c("Other", "RBC")) +
  scale_alpha_manual(values=c(0.3, 1), labels=c("Other", "RBC")) +
  guides(alpha="none") +
  ggtitle("APR+PCA+UMAP")

p <- ggarrange(p_apr_2, p_apr_umap_2, nrow=1, ncol=2,common.legend=TRUE)

ggsave(p, filename="../plots/blish_apr_example.png",
       width=10.8, height=7.24, units="in")


load("../data/simdata_time.RData")
load("../data/simdata_LL.RData")


library(tidyverse)
## Plotting function for loglike vs runtime comparison


#cutoff <- min(which(abs(rel.diff) < 10^{-4}))

df.1 <- data.frame(time=Time[[1]],
                   Method="scGBM-full",
                   ll=LL[[1]])


df.2 <- data.frame(time=Time[[2]],
                   Method="scGBM-proj",
                   ll=LL[[2]])


df.3 <- data.frame(time=Time[[3]],
                   Method="GLM-PCA (AvaGrad)",
                   ll=LL[[3]])

df.4 <- data.frame(time=Time[[4]],
                   Method="GLM-PCA (Fisher)",
                   ll=LL[[4]])

df.5 <- data.frame(time=Time[[5]],
                   Method="GLM-PCA (SGD)",
                   ll=LL[[5]])



df <- rbind(df.1,df.2,df.3,df.4,df.5) %>% as.data.frame


p <- ggplot(data=df,aes(x=time/3600,y=ll,color=Method))
p <- p + geom_point() + geom_line()
#p <- p + geom_segment(x=log10(Time[[2]][1]/3600),xend=100,y=LL[[2]][1],yend=LL[[2]][1],
#color=hue_pal()(5)[5],linetype="dashed")
p <- p + scale_x_log10(
  breaks = scales::trans_breaks("log10", function(x) 10^x),
  labels = scales::trans_format("log10", scales::math_format(10^.x)),
  limits=c(10^{1.2},10^{3.5})/3600
)
p <- p + xlab("") + ylab("")+
  ylim(2.24*10^8, NA)
p <- p + theme_bw()
p <- p + annotation_logticks(sides = 'b')
p <- p + ggtitle("LLV Sim (J=100,000; I = 1000)")
p <- p + guides(color="none")
p_simdata <- p



load("../data/simdatalv5_time.RData")
load("../data/simdatalv5_LL.RData")


library(tidyverse)
## Plotting function for loglike vs runtime comparison


#cutoff <- min(which(abs(rel.diff) < 10^{-4}))

df.1 <- data.frame(time=Time[[1]],
                   Method="scGBM-full",
                   ll=LL[[1]])


df.2 <- data.frame(time=Time[[2]],
                   Method="scGBM-proj",
                   ll=LL[[2]])


df.3 <- data.frame(time=Time[[3]],
                   Method="GLM-PCA (AvaGrad)",
                   ll=LL[[3]])

df.4 <- data.frame(time=Time[[4]],
                   Method="GLM-PCA (Fisher)",
                   ll=LL[[4]])

df.5 <- data.frame(time=Time[[5]],
                   Method="GLM-PCA (SGD)",
                   ll=LL[[5]])



df <- rbind(df.1,df.2,df.3,df.4,df.5) %>% as.data.frame


p <- ggplot(data=df,aes(x=time/3600,y=ll,color=Method))
p <- p + geom_point() + geom_line()
#p <- p + geom_segment(x=log10(Time[[2]][1]/3600),xend=100,y=LL[[2]][1],yend=LL[[2]][1],
#color=hue_pal()(5)[5],linetype="dashed")
p <- p + scale_x_log10(
  breaks = scales::trans_breaks("log10", function(x) 10^x),
  labels = scales::trans_format("log10", scales::math_format(10^.x)),
  limits=c(10^{1.2},10^{3.5})/3600
)
p <- p + xlab("") + ylab("")+
  ylim(2.5*10^8, NA)
p <- p + theme_bw()
p <- p + annotation_logticks(sides = 'b')
p <- p + ggtitle("HLV Simulation (J=100,000; I = 1000)")
p <- p + guides(color="none")
p_simdatalv5 <- p


library(ggpubr)

p_ll <- ggarrange(pblish, p_tenximmune,
               p_simdata, p_simdatalv5, nrow=2, ncol=2, common.legend = TRUE,
               legend="bottom")

ggsave(p_ll, filename="../plots/all_runtime_plots.png")



### Plotting
res <- readRDS("../../simulation_accuracy/data/factor_correlation.RDS")
res2 <- readRDS("../../simulation_accuracy/data/mse.RDS")

df <- reshape2::melt(res)


I <- 10^3
J <- 10^4
M <- 10

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


#res2 <- sqrt(I*J)*res2
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

p <- ggarrange(p, p_ll, nrow=2, heights=c(1.5,2), labels=c("a","b"))

ggsave(p, filename="../plots/glmpca_comparison.png",
       width=11.1, height=11.7)














load("../data/zhengmix_Time_gbmInit.RData")
load("../data/zhengmix_LL_gbmInit.RData")


library(tidyverse)
## Plotting function for loglike vs runtime comparison


#cutoff <- min(which(abs(rel.diff) < 10^{-4}))


df.3 <- data.frame(time=Time[[1]],
                   Method="GLM-PCA (AvaGrad)",
                   ll=LL[[1]])

df.4 <- data.frame(time=Time[[2]],
                   Method="GLM-PCA (Fisher)",
                   ll=LL[[2]])

df.5 <- data.frame(time=Time[[3]],
                   Method="GLM-PCA (SGD)",
                   ll=LL[[3]])


load("../data/zhengmix_Time.RData")
load("../data/zhengmix_LL.RData")

df.1 <- data.frame(time=Time[[1]],
                   Method="scGBM-full",
                   ll=LL[[1]])


df.2 <- data.frame(time=rowMeans(Time[[2]]),
                   Method="scGBM-proj",
                   ll=rowMeans(LL[[2]]))



df <- rbind(df.1, df.2, df.3,df.4,df.5) %>% as.data.frame


p <- ggplot(data=df,aes(x=time/3600,y=ll,color=Method))
p <- p + geom_point() + geom_line()
#p <- p + geom_segment(x=log10(Time[[2]][1]/3600),xend=100,y=LL[[2]][1],yend=LL[[2]][1],
#color=hue_pal()(5)[5],linetype="dashed")
p <- p + scale_x_log10(
  breaks = scales::trans_breaks("log10", function(x) 10^x),
  labels = scales::trans_format("log10", scales::math_format(10^.x)),
  #limits=c(10,10^{3})/3600
)
p <- p + xlab("") + ylab("")
p <- p + theme_bw()
p <- p + annotation_logticks(sides = 'b')
p <- p + ggtitle("10X Immune (J=3,994; I=6,049)")
#p <- p + guides(color="none")
p  <- p
p_tenximmune <- p

ggsave(p, filename="../plots/10x_with_gbm_init.png")



