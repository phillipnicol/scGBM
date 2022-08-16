plot.gbm <- function(gbm,dims=c(1,2),cluster=NULL,
                     se_V=NULL,return.gg=FALSE) {
  V <- out$V
  if(is.null(se_V)) {
    df <- data.frame(x=V[,dims[1]],y=V[,dims[2]])
    if(is.null(cluster)) {
      df$cluster <- rep(1,nrow(df))
    } else{
      df$cluster <- cluster
    }
    p <- ggplot(data=df,aes(x=x,y=y,color=cluster))
    p <- p + geom_point()
    print(p)
  } else {
    df <- data.frame(x0=V[,dims[1]],y0=V[,dims[2]])
    if(is.null(cluster)) {
      df$cluster <- rep(1,nrow(df))
    } else{
      df$cluster <- cluster
    }
    df$a <- se_V[,dims[1]]; df$b <- se_V[,dims[2]]
    df$angle <- rep(0,nrow(df))
    p <- ggplot(data=df,aes(x0=x0,y0=y0,fill=cluster,
                            a=a,b=b,angle=angle))
    p <- p + geom_ellipse(alpha=0.5)
    print(p)
  }
  if(return.gg) {
    return(p)
  }
}

