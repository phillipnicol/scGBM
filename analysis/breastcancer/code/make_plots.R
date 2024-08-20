
library(tidyverse)

meta <- read.csv("../../data/breastcancer/metadata.csv")

V <- readRDS("../data/V.RDS")

p <- data.frame(x=V[,1], y=V[,2], color=meta$celltype_major) |>
  ggplot(aes(x=x,y=y,color=color)) +
  geom_point(size=0.25) +
  theme_bw()

H.table <- readRDS("../data/H.table.minor.RDS") |> as.data.frame()




p <- pheatmap::pheatmap(H.table,legend=TRUE, color=colorRampPalette(c("white","red"))(100),
                   breaks=seq(0,1,by=0.01),
                   rownames=TRUE,
                   colnames=TRUE,
                   cluster_rows = ifelse(nrow(H.table) > 1, TRUE,FALSE),
                   cluster_cols = ifelse(nrow(H.table) > 1, TRUE,FALSE))

Vu <- readRDS("../data/gbm.umap.RDS")

p <- data.frame(x=Vu[,1], y=Vu[,2], color=ifelse(meta$celltype_minor %in% c("B cells Memory",
                                                                            "B cells Naive",
                                                                            "NK cells",
                                                                            "NKT cells",
                                                                            "Cycling T-cells",
                                                                            "T cells CD8+",
                                                                            "T cells CD4+"),
                                                 meta$celltype_minor,
                                                 "")) |>
  ggplot(aes(x=x,y=y,color=color)) +
  geom_point(size=0.25) +
  theme_bw()



p <- data.frame(x=Vu[,1], y=Vu[,2], color=ifelse(meta$celltype_minor %in% c("Luminal Progenitors",
                                                                            "Myoepithelial",
                                                                            "Plasmablasts",
                                                                            "Cycling PVL"),
                                                 meta$celltype_minor,
                                                 "")) |>
  ggplot(aes(x=x,y=y,color=color)) +
  geom_point(size=0.25) +
  theme_bw()

colors <- rainbow(12)

p_high_cci_low_icc <- data.frame(x=Vu[,1], y=Vu[,2], color=ifelse(meta$celltype_minor %in% c("Luminal Progenitors",
                                                                            "Myoepithelial",
                                                                            "Plasmablasts",
                                                                            "Cycling PVL"),
                                                 meta$celltype_minor,
                                                 " ")) |>
  ggplot(aes(x=x,y=y,color=color, alpha=color)) +
  scale_color_manual(values = c(" " = "grey",
                                "Luminal Progenitors" = colors[1],
                                "Myoepithelial" = colors[2],
                                "Plasmablasts" = colors[3],
                                "Cycling PVL" = colors[4]),
                     breaks=c("Luminal Progenitors",
                              "Myoepithelial",
                              "Plasmablasts",
                              "Cycling PVL")) +
  scale_alpha_manual(values = c(" " = 0.3,
                                "Luminal Progenitors" = 1,
                                "Myoepithelial" = 1,
                                "Plasmablasts" = 1,
                                "Cycling PVL" = 1)) +
  guides(alpha="none") + labs(color="Cell type") +
  geom_point(size=0.25) + ggtitle("High CCI + Low inter-cluster CCI") +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2")+ theme(legend.position = "bottom")

p_high_cci_high_icc <- data.frame(x=Vu[,1], y=Vu[,2], color=ifelse(meta$celltype_minor %in% c("Endothelial Lymphatic LYVE1",
                                                                            "Endothelial ACKR1",
                                                                            "Endothelial RGS5",
                                                                            "Endothelial CXCL12"),
                                                 meta$celltype_minor,
                                                 " ")) |>
  ggplot(aes(x=x,y=y,color=color, alpha=color)) +
  scale_color_manual(values = c(" " = "grey",
                                "Endothelial Lymphatic LYVE1" = colors[5],
                                "Endothelial ACKR1" = colors[6],
                                "Endothelial RGS5" = colors[7],
                                "Endothelial CXCL12" = colors[8]),
                     breaks=c("Endothelial Lymphatic LYVE1",
                              "Endothelial ACKR1",
                              "Endothelial RGS5",
                              "Endothelial CXCL12")) +
  scale_alpha_manual(values = c(" " = 0.3,
                                "Endothelial Lymphatic LYVE1" = 1,
                                "Endothelial ACKR1" = 1,
                                "Endothelial RGS5" = 1,
                                "Endothelial CXCL12" = 1)) +
  guides(alpha="none") + labs(color="Cell type") +
  geom_point(size=0.25) + ggtitle("High CCI + High inter-cluster CCI") +
  xlab("UMAP1") + ylab("UMAP2") + theme_bw()+ theme(legend.position = "bottom")


p_low_cci <- data.frame(x=Vu[,1], y=Vu[,2], color=ifelse(meta$celltype_minor %in% c("Mature Luminal",
                                                                            "Cancer Basal SC",
                                                                            "Cancer Cycling",
                                                                            "Cancer Her2 SC"),
                                                 meta$celltype_minor,
                                                 " ")) |>
  ggplot(aes(x=x,y=y,color=color,alpha=color)) +
  geom_point(size=0.25) +
  scale_color_manual(values = c(" " = "grey",
                                "Mature Luminal" = colors[9],
                                "Cancer Basal SC" = colors[10],
                                "Cancer Cycling" = colors[11],
                                "Cancer Her2 SC" = colors[12]),
                     breaks = c("Mature Luminal",
                                "Cancer Basal SC",
                                "Cancer Cycling",
                                "Cancer Her2 SC")) +
  scale_alpha_manual(values = c(" " = 0.3,
                                "Mature Luminal" = 1,
                                "Cancer Basal SC" = 1,
                                "Cancer Cycling" = 1,
                                "Cancer Her2 SC" = 1)) +
  xlab("UMAP1") + ylab("UMAP2") +
  ggtitle("Low CCI") +
  guides(alpha="none") + labs(color="Cell type") +
  theme_bw()

p <- ggarrange(p_high_cci_low_icc, p_high_cci_high_icc, p_low_cci, ncol=3)


cci <- readRDS("../data/cci_minor.RDS")
