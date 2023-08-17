loadings.volcano <- function(gbm, dim=1) {
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
}

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rownames(out$loadings)

ensembl = useMart(biomart = "ensembl", dataset="hsapiens_gene_ensembl")
res <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),
      filters="ensembl_gene_id",
      values=genes,
      mart=ensembl)

for(i in 1:nrow(out$loadings)) {
  ix <- which(res$ensembl_gene_id == rownames(out$loadings)[i])
  if(length(ix) > 0) {
    rownames(out$loadings)[i] <- res$hgnc_symbol[i]
  }
}
