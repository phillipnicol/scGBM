#' @export
#'
#' @title Plot the factor loadings
#'
#' @description A volcano plot shows which genes are driving
#' a particular scGBM factor
#'
#' @param gbm A list that is the return value of \code{gbm.sc}
#' @param dim Which latent factor to plot
#' @param return.plot Should ggplot object containing the plot be returned?
#'
#' @return If `return.plot` is TRUE, a `ggplot2` object is returned.
#'
#' @author Phillip B. Nicol <philnicol740@gmail.com>
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
