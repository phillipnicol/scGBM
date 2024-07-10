

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

p <- ggarrange(pblish, p_tenximmune,
               p_simdata, p_simdatalv5, nrow=2, ncol=2, common.legend = TRUE)

ggsave(p, filename="../plots/all_runtime_plots.png")

p_isba <- ggarrange(pblish, p_tenximmune,
               p_simdata, p_simdatalv5, nrow=1, ncol=4, common.legend = TRUE)
p_isba <- annotate_figure(p_isba,
                          left = text_grob("Out-of-sample LL", rot = 90,size=20),
                          bottom = text_grob("Wall time (hours)",size=20))


ggsave(p_isba, filename="../plots/all_runtime_plots_isba.pdf",
       width=20.7, height=5.96, units="in")
