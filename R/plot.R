#' @import ggplot2
NULL

#' @export
#'
#' @title Plot the results of scGBM
#'
#' @param gbm An object which is the output of \code{\link{gbm.sc}}
#' @param dims Which dimensions (in the latent factors) to plot
#' @param se If TRUE, plot the uncertainty estimates around each point
#' @param return.gg If TRUE, return the ggplot object
#' @param plot If TRUE, plot the \code{ggplot} object
#'
#' @return If \code{return.gg} is TRUE, then a \code{ggplot} object is
#' returned.
#'
#' @author Phillip B. Nicol <philnicol740@gmail.com>
plot_gbm <- function(gbm,dims=c(1,2),cluster=NULL,
                     se=FALSE,return.gg=FALSE, plot=TRUE) {
  V <- out$scores
  if(!se) {
    df <- data.frame(x=V[,dims[1]],y=V[,dims[2]])
    if(is.null(cluster)) {
      df$cluster <- rep(1,nrow(df))
    } else{
      df$cluster <- cluster
    }
    p <- ggplot(data=df,aes(x=x,y=y,color=cluster))
    p <- p + geom_point()
    p <- p + xlab(paste0("GBM-",dims[1])) + ylab(paste0("GBM-",dims[2]))
    print(p)
  } else {
    se_V <- out$se_scores
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
    p <- p + xlab(paste0("GBM-",dims[1])) + ylab(paste0("GBM-",dims[2]))
    p <- p + ggforce::geom_ellipse(alpha=0.5)
    print(p)
  }
  if(return.gg) {
    return(p)
  }
}

