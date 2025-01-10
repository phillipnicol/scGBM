
library(tidyverse)
library(ggpubr)

meta <- read.csv("../../data/breastcancer/metadata.csv")

V <- readRDS("../data/V.RDS")

p.embedding <- data.frame(x=V[,1], y=V[,2], color=meta$celltype_major) |>
  ggplot(aes(x=x,y=y,color=color)) +
  geom_point(size=0.25) +
  theme_bw() +
  xlab("GBM1") + ylab("GBM2")

H.table <- readRDS("../data/H.table.minor.RDS") |> as.data.frame()


colors <- rainbow(12)

p.heatmap <- pheatmap::pheatmap(H.table,legend=TRUE, color=colorRampPalette(c("white","red"))(100),
                   breaks=seq(0,1,by=0.01),
                   rownames=TRUE,
                   colnames=TRUE,
                   cluster_rows = ifelse(nrow(H.table) > 1, TRUE,FALSE),
                   cluster_cols = ifelse(nrow(H.table) > 1, TRUE,FALSE))

library(grid)
my.colors <- rep("black", unique(meta$celltype_minor) |> length())
names(my.colors) <- unique(meta$celltype_minor)[p.heatmap$tree_row$order]


my.colors["Luminal Progenitors"] <- colors[1]
my.colors["Myoepithelial"] <- colors[4]
my.colors["Plasmablasts"] <- colors[7]
my.colors["Cycling PVL"] <- colors[10]

my.colors["Endothelial Lymphatic LYVE1"] <- colors[2]
my.colors["Endothelial ACKR1"] <- colors[5]
my.colors["Endothelial RGS5"] <- colors[8]
my.colors["Endothelial CXCL12"] <- colors[11]

my.colors["Mature Luminal"] <- colors[3]
my.colors["Cancer Basal SC"] <- colors[6]
my.colors["Cancer Cycling"] <- colors[9]
my.colors["Cancer Her2 SC"] <- colors[12]

p.heatmap$gtable$grobs[[4]]$gp=gpar(col=my.colors)
p.heatmap$gtable$grobs[[5]]$gp=gpar(col=my.colors)


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


p_high_cci_low_icc <- data.frame(x=Vu[,1], y=Vu[,2], color=ifelse(meta$celltype_minor %in% c("Luminal Progenitors",
                                                                            "Myoepithelial",
                                                                            "Plasmablasts",
                                                                            "Cycling PVL"),
                                                 meta$celltype_minor,
                                                 " ")) |>
  ggplot(aes(x=x,y=y,color=color, alpha=color)) +
  scale_color_manual(values = c(" " = "grey",
                                "Luminal Progenitors" = colors[1],
                                "Myoepithelial" = colors[4],
                                "Plasmablasts" = colors[7],
                                "Cycling PVL" = colors[10]),
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
  geom_point(size=0.25) + ggtitle("High CCI + Low inter-CCI") +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2")+ theme(legend.position = "bottom") +
  guides(color="none")

p_high_cci_high_icc <- data.frame(x=Vu[,1], y=Vu[,2], color=ifelse(meta$celltype_minor %in% c("Endothelial Lymphatic LYVE1",
                                                                            "Endothelial ACKR1",
                                                                            "Endothelial RGS5",
                                                                            "Endothelial CXCL12"),
                                                 meta$celltype_minor,
                                                 " ")) |>
  ggplot(aes(x=x,y=y,color=color, alpha=color)) +
  scale_color_manual(values = c(" " = "grey",
                                "Endothelial Lymphatic LYVE1" = colors[2],
                                "Endothelial ACKR1" = colors[5],
                                "Endothelial RGS5" = colors[8],
                                "Endothelial CXCL12" = colors[11]),
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
  geom_point(size=0.25) + ggtitle("High CCI + High inter-CCI") +
  xlab("UMAP1") + ylab("UMAP2") + theme_bw()+ theme(legend.position = "bottom") +
  guides(color="none")


p_low_cci <- data.frame(x=Vu[,1], y=Vu[,2], color=ifelse(meta$celltype_minor %in% c("Mature Luminal",
                                                                            "Cancer Basal SC",
                                                                            "Cancer Cycling",
                                                                            "Cancer Her2 SC"),
                                                 meta$celltype_minor,
                                                 " ")) |>
  ggplot(aes(x=x,y=y,color=color,alpha=color)) +
  geom_point(size=0.25) +
  scale_color_manual(values = c(" " = "grey",
                                "Mature Luminal" = colors[3],
                                "Cancer Basal SC" = colors[6],
                                "Cancer Cycling" = colors[9],
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
  theme_bw() + guides(color="none")

p <- ggarrange(p_high_cci_low_icc, p_high_cci_high_icc, p_low_cci, ncol=3)


p.full <- ggarrange(p.embedding,
          p.heatmap[[4]],
          p,
          labels=c("a","b","c"),
          nrow=3,ncol=1,
          heights=c(1,2,1))

ggsave(p.full, filename="../plots/gbm_umap_cci_breastcancer.png",
       width=10.3, height=16, units="in")

#Editable for illustrator
#ggsave(p.full, filename="../plots/gbm_umap_cci_breastcancer.pdf",
#       width=10.3, height=16, units="in")


## Look into Cycling PVL

p.list <- list()
for(i in 1:5) {
  p <- data.frame(x=V[,2*i-1], y=V[,2*i], color=ifelse(meta$celltype_minor == "Cycling PVL", "Cycling PVL", "Other")) |>
    arrange(color=="Cycling PVL") |>
    ggplot(aes(x=x,y=y,color=color)) +
    geom_point(size=0.5) +
    scale_color_manual(labels=c("Cycling PVL", "Other"), values=c("firebrick", "grey90")) +
    theme_bw() +
    labs(color="") +
    xlab(paste0("GBM", 2*i - 1)) + ylab(paste0("GBM", 2*i))

  p.list[[i]] <- p
}

p <- ggarrange(p.list[[1]],
          p.list[[2]],
          p.list[[3]],
          p.list[[4]],
          p.list[[5]],
          common.legend = TRUE,
          nrow=1)

ggsave(p, filename="../plots/cycling_pvl_gbm.png",
       width=13, height=3.92, units="in")

Vu <- readRDS("../data/aprUMAP.RDS")



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
  geom_point(size=0.25) + ggtitle("High CCI + Low inter-CCI") +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2")+ theme(legend.position = "bottom") +
  guides(color="none")


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
  geom_point(size=0.25) + ggtitle("High CCI + High inter-CCI") +
  xlab("UMAP1") + ylab("UMAP2") + theme_bw()+ theme(legend.position = "bottom") +
  guides(color="none")


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
  theme_bw() + guides(color="none")

p <- ggarrange(p_high_cci_low_icc, p_high_cci_high_icc, p_low_cci, ncol=3)


##LOG+PCA

Vu <- readRDS("../data/pca.umap.RDS")



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
  geom_point(size=0.25) + ggtitle("High CCI + Low inter-CCI") +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2")+ theme(legend.position = "bottom") +
  guides(color="none")

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
  geom_point(size=0.25) + ggtitle("High CCI + High inter-CCI") +
  xlab("UMAP1") + ylab("UMAP2") + theme_bw()+ theme(legend.position = "bottom") +
  guides(color="none")


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
  theme_bw() + guides(color="none")

p.log.pca <- ggarrange(p_high_cci_low_icc, p_high_cci_high_icc, p_low_cci, ncol=3)
p.log.pca <- annotate_figure(p.log.pca, top = text_grob("Log+PCA+UMAP",
                                                            color = "black", face = "bold", size = 14))



Vu <- readRDS("../data/aprUMAP.RDS")



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
  geom_point(size=0.25) + ggtitle("High CCI + Low inter-CCI") +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2")+ theme(legend.position = "bottom") +
  guides(color="none")


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
  geom_point(size=0.25) + ggtitle("High CCI + High inter-CCI") +
  xlab("UMAP1") + ylab("UMAP2") + theme_bw()+ theme(legend.position = "bottom") +
  guides(color="none")


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
  theme_bw() + guides(color="none")

p.apr <- ggarrange(p_high_cci_low_icc, p_high_cci_high_icc, p_low_cci, ncol=3)
p.apr <- annotate_figure(p.apr, top = text_grob("APR+PCA+UMAP",
                                                            color = "black", face = "bold", size = 14))



#LOG+SCALE+PCA


Vu <- readRDS("../data/log.scale.pca.UMAP.RDS")


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
  geom_point(size=0.25) + ggtitle("High CCI + Low inter-CCI") +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2")+ theme(legend.position = "bottom") +
  guides(color="none")


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
  geom_point(size=0.25) + ggtitle("High CCI + High inter-CCI") +
  xlab("UMAP1") + ylab("UMAP2") + theme_bw()+ theme(legend.position = "bottom") +
  guides(color="none")


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
  theme_bw() + guides(color="none")

p.scale.pca <- ggarrange(p_high_cci_low_icc, p_high_cci_high_icc, p_low_cci, ncol=3)
p.scale.pca <- annotate_figure(p.scale.pca, top = text_grob("Log+Scale+PCA+UMAP",
                                      color = "black", face = "bold", size = 14))


p.default.embed <- ggarrange(p.scale.pca,
                             p.log.pca,
                             p.apr,
                             nrow=3,ncol=1)

ggsave(p.default.embed, filename="../plots/breastcancer_default_embedding.png",
       width=2*7.89, height=2*5.99, units="in")

