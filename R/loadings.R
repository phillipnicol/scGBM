
plot.loadings <- function(fit, dim=1, n=5) {
  Um <- fit$U[,dim]

  topn <- order(Um, decreasing=TRUE)[1:n]
  bottomn <- order(Um)[1:n]

  if(is.null(rownames(fit$U))) {
    rownames(fit$U) <- 1:nrow(fit$U)
  }

  df <- data.frame(y = as.character(rownames(fit$U)[c(topn,bottomn)]),
                   x = Um[c(topn,bottomn)])

  p <- ggplot(data=df,aes(x=x,y=y))
  p <- p + geom_point()
}
