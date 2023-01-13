


fn <- function(V,nstart=25,centers=5) {
  kmeans(V,nstart=nstart,centers=centers)$cluster
}

CCI <- function(V,
                se_V,
                cluster.orig,
                reps=100,
                cluster.fn=fn,
                ...) {

  nc <- length(table(cluster.orig))
  cci <- matrix(0,nrow=reps,ncol=nc)

  J <- nrow(V)
  for(i in 1:reps) {
    print(i)
    M <- length(out$D)
    e <- matrix(rnorm(n=J*M,sd=se_V),nrow=J,ncol=M)

    V.new <- V+e
    cluster <- fn(V.new,...)

    for(k in 1:nc) {
      ixs <- which(cluster.orig==unique(cluster.orig)[k])
      new.ixs <- cluster[ixs]
      cci[i,k] <- sum(choose(table(new.ixs),2))/choose(length(ixs),2)
    }
  }

  ##Plotting
  colnames(cci) <- unique(cluster.orig)
  df <- reshape2::melt(cci)
  df$Var2 <- as.factor(df$Var2)

  p <- ggplot(data=df,aes(x=Var2,y=value,fill=Var2))
  p <- p + geom_boxplot(alpha=0.5)
  p <- p + theme_bw()
  p <- p + xlab("Cluster") + ylab("CCI")
  p <- p + ggtitle("Posterior of Cluster Confidence Index")
  p  <- p + guides(fill="none")
  p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  p

  return(p)
}


CCI.heatmap <- function(V,
                        se_V,
                        cluster.orig,
                        reps=100,
                        cluster.fn=fn,
                        ...) {

  nc <- length(table(cluster.orig))
  cci <- array(0,dim=c(nc,nc,reps))
  H.table <- matrix(0,nrow=nc,ncol=nc)
  J <- nrow(V)
  for(i in 1:reps) {
    print(i)
    M <- length(out$D)
    e <- matrix(rnorm(n=J*M,sd=se_V),nrow=J,ncol=M)

    V.new <- V+e
    cluster <- fn(V.new,...)

    for(k in 1:nc) {
      for(j in 1:nc) {
        kxs <- which(cluster.orig==unique(cluster.orig)[k])
        jxs <- which(cluster.orig==unique(cluster.orig)[j])

        new.jxs <- cluster[jxs]
        new.kxs <- cluster[kxs]

        new.jxs <- factor(new.jxs,levels=c(1:nc))
        new.kxs <- factor(new.kxs,levels=c(1:nc))

        tab.jxs <- table(new.jxs)
        tab.kxs <- table(new.kxs)

        overlap <- sum(tab.jxs*tab.kxs)
        max.overlap <- length(kxs)*length(jxs)

        cci[k,j,i] <- overlap/max.overlap

      }
    }
  }

  for(i in 1:nc) {
    for(j in 1:nc) {
      H.table[i,j] <- mean(cci[i,j,])
    }
  }

  H.table <- data.frame(H.table)

  x1 <- unique(cluster.orig)
  x2 <- unique(cluster.orig)

  my.grid <- expand.grid(x1,x2)
  v <- c()

  cntr <- 1
  for(i in 1:nc) {
    for(j in 1:nc) {
      v <- c(v,H.table[i,j])
    }
  }

  df <- cbind(my.grid,v)
  colnames(df) <- c("x","y","z")

  p <- ggplot(data=df,aes(x=x,y=y,fill=z)) + geom_tile()
  p <- p + scale_fill_gradient(low="blue", high="red",
                                limits = c(0,1))
  p <- p + theme_bw()
  p <- p + xlab("Cluster") + ylab("Cluster")
  p <- p + labs(fill="Posterior mean")

  return(p)
}

