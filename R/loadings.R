loadings.volcano <- function(gbm, dim=1, return.plot=FALSE) {
  Um <- gbm$loadings[,dim]

  wald <- (Um/gbm$se_loadings[,dim])^2
  pval <- pchisq(wald, df=1, lower.tail=FALSE)
  pval[pval < 10^{-15}] <- 10^{-15}

  I <- nrow(gbm$loadings)

  xmax <- max(abs(Um))

  res <- data.frame(log2FC=Um,pvalue=pval)

  p <- EnhancedVolcano::EnhancedVolcano(res,x="log2FC",y="pvalue",
                  lab=rownames(out$loadings),FCcutoff=2/(sqrt(I)),
                  xlim=c(-1.05*xmax, 1.05*xmax),
                  xlab="Loading",legendPosition="none",
                  ylim=c(0,15),
                  title=NULL,
                  subtitle=paste0("Factor ", dim))
  p <- p + geom_vline(xintercept=sqrt(2*log(2*I)/I), color="red",
                      linetype="dashed")
  p <- p + geom_vline(xintercept=-sqrt(2*log(2*I)/I),color="red",
                      linetype="dashed")
  p

  if(return.plot) {
    return(p)
  }
}
