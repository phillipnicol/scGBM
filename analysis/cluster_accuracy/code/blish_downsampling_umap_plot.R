setwd(here::here("analysis/cluster_accuracy/code"))


apr <- readRDS("../data/covid19_subsample_apr_umap.RDS")

glmpca <- readRDS("../data/covid19_subsample_glmpca_umap.RDS")

lpca_scale <- readRDS("../data/covid19_subsample_lpca_scale_umap.RDS")

lpca <- readRDS("../data/covid19_subsample_lpca_umap.RDS")

sct <- readRDS("../data/covid19_subsample_sct_scale_umap.RDS")

gbm <- readRDS("../data/covid19_subsample_scGBM_umap.RDS")


meta <- readRDS("../../data/blish/blish_meta.RDS")


size <- 0.25

library(tidyverse)

p.apr <- apr |> as.data.frame() |>
  mutate(cell_type = meta$cell.type.coarse) |>
  ggplot(aes(x=umap_1, y= umap_2, color=cell_type)) +
  geom_point(size=size) + theme_bw() +
  xlab("UMAP 1") + ylab("UMAP 2") +
  labs(color= "Cell type") +
  ggtitle("APR+PCA+UMAP")

p.glmpca <- glmpca |> as.data.frame() |>
  mutate(cell_type = meta$cell.type.coarse) |>
  ggplot(aes(x=umap_1, y= umap_2, color=cell_type)) +
  geom_point(size=size) + theme_bw() +
  xlab("UMAP 1") + ylab("UMAP 2") +
  labs(color= "Cell type") +
  ggtitle("GLMPCA+UMAP")

p.lpcascale <- lpca_scale |> as.data.frame() |>
  mutate(cell_type = meta$cell.type.coarse) |>
  ggplot(aes(x=umap_1, y= umap_2, color=cell_type)) +
  geom_point(size=size) + theme_bw() +
  xlab("UMAP 1") + ylab("UMAP 2") +
  labs(color= "Cell type") +
  ggtitle("Log+Scale+PCA+UMAP")


p.lpca <- lpca |> as.data.frame() |>
  mutate(cell_type = meta$cell.type.coarse) |>
  ggplot(aes(x=umap_1, y= umap_2, color=cell_type)) +
  geom_point(size=size) + theme_bw() +
  xlab("UMAP 1") + ylab("UMAP 2") +
  labs(color= "Cell type") +
  ggtitle("LOG+PCA+UMAP")


p.sct <- sct |> as.data.frame() |>
  mutate(cell_type = meta$cell.type.coarse) |>
  ggplot(aes(x=umap_1, y= umap_2, color=cell_type)) +
  geom_point(size=size) + theme_bw() +
  xlab("UMAP 1") + ylab("UMAP 2") +
  labs(color= "Cell type") +
  ggtitle("SCT+PCA+UMAP")


p.gbm <- gbm |> as.data.frame() |>
  mutate(cell_type = meta$cell.type.coarse) |>
  ggplot(aes(x=umap_1, y= umap_2, color=cell_type)) +
  geom_point(size=size) + theme_bw() +
  xlab("UMAP 1") + ylab("UMAP 2") +
  labs(color= "Cell type") +
  ggtitle("scGBM+UMAP")




library(ggpubr)

p <- ggarrange(p.gbm, p.apr,
          p.sct, p.lpcascale,
          p.lpca, p.glmpca,
          nrow=3, ncol=2,
          common.legend = TRUE,
          legend="top")



ggsave(p, filename="../plots/blish_downsampled_umap.png", width=15, height=10, units="in")


