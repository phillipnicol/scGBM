
#Zhengmix
load("../data/zhengmix_LL.Rdata")
load("../data/zhengmix_time.Rdata")


library(tidyverse)
## Plotting function for loglike vs runtime comparison

ixs <- which(!is.infinite(LL[[1]]))
df.1 <- data.frame(time=Time[[1]][ixs],
                   method="scGBM-full",
                   ll=LL[[1]][ixs])


df.2 <- data.frame(time=Time[[2]],
                   method="scGBM-proj",
                   ll=LL[[2]])


df.3 <- data.frame(time=Time[[3]],
                   method="GLM-PCA (AvaGrad)",
                   ll=LL[[3]])

df.4 <- data.frame(time=Time[[4]],
                   method="GLM-PCA (Fisher)",
                   ll=LL[[4]])

df.5 <- data.frame(time=Time[[5]],
                   method="GLM-PCA (SGD)",
                   ll=LL[[5]])



df <- rbind(df.1,df.2,df.3,df.4,df.5) %>% as.data.frame


p <- ggplot(data=df,aes(x=time/3600,y=ll,color=method))
p <- p + geom_point() + geom_line()
#p <- p + geom_segment(x=log10(Time[[2]][1]/3600),xend=100,y=LL[[2]][1],yend=LL[[2]][1],
#color=hue_pal()(5)[5],linetype="dashed")
p <- p + scale_x_log10(
  breaks = scales::trans_breaks("log10", function(x) 10^x),
  labels = scales::trans_format("log10", scales::math_format(10^.x)),
  limits=c(1,10^{3})/3600
)
p <- p + xlab("") + ylab("")
p <- p + theme_bw()
p <- p + annotation_logticks(sides = 'b')
p <- p + ggtitle("10X Immune (J=3,994)")
p <- p + guides(color="none")
p_tenximmune <- p




