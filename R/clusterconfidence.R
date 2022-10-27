


fn <- function(V) {
  kmeans(V,centers=K,nstart=25)$cluster
}

cluster.confidence <- function(V,
                               se_V,
                               cluster.orig,
                               reps=100,
                               null.dist=FALSE) {

  nc <- length(table(cluster.orig))
  cci <- matrix(0,nrow=reps,ncol=nc)
  cci.null <- matrix(0,nrow=reps,ncol=nc)
  for(i in 1:reps) {
    print(i)
    M <- length(out$D)
    e <- matrix(rnorm(n=J*M,sd=se_V),nrow=J,ncol=M)

    V.new <- V+e
    cluster <- kmeans(V.new,centers=nc,nstart=25)

    for(k in 1:nc) {
      ixs <- which(cluster.orig==unique(cluster.orig)[k])
      new.ixs <- cluster$cluster[ixs]
      cci[i,k] <- sum(choose(table(new.ixs),2))/choose(length(ixs),2)
    }

    if(null.dist) {
      V.null <- e
      cluster <- kmeans(V.new,centers=nc,nstart=25)

      for(k in 1:nc) {
        ixs <- which(cluster.orig==unique(cluster.orig)[k])
        new.ixs <- cluster$cluster[ixs]
        cci.null[i,k] <- sum(choose(table(new.ixs),2))/choose(length(ixs),2)
      }
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
}


cluster.confidence2 <- function(V,
                               se_V,
                               cluster.orig,
                               reps=100) {

  nc <- length(table(cluster.orig))
  cci <- array(0,dim=c(nc,nc,reps))
  cci.null <- cci
  for(i in 1:reps) {
    print(i)
    M <- length(out$D)
    e <- matrix(rnorm(n=J*M,sd=se_V),nrow=J,ncol=M)

    V.new <- V+e
    cluster <- kmeans(V.new,centers=nc,nstart=25)

    for(k in 1:nc) {
      for(j in 1:nc) {
        kxs <- which(cluster.orig==unique(cluster.orig)[k])
        jxs <- which(cluster.orig==unique(cluster.orig)[j])

        new.jxs <- cluster$cluster[jxs]
        new.kxs <- cluster$cluster[kxs]

        new.jxs <- factor(new.jxs,levels=c(1:nc))
        new.kxs <- factor(new.kxs,levels=c(1:nc))

        tab.jxs <- table(new.jxs)
        tab.kxs <- table(new.kxs)

        overlap <- sum(tab.jxs*tab.kxs)
        max.overlap <- length(kxs)*length(jxs)

        cci[k,j,i] <- overlap/max.overlap

      }
    }

    V.null <- e
    cluster <- kmeans(V.null,centers=nc,nstart=25)

    for(k in 1:nc) {
      for(j in 1:nc) {
        kxs <- which(cluster.orig==unique(cluster.orig)[k])
        jxs <- which(cluster.orig==unique(cluster.orig)[j])

        new.jxs <- cluster$cluster[jxs]
        new.kxs <- cluster$cluster[kxs]

        new.jxs <- factor(new.jxs,levels=c(1:nc))
        new.kxs <- factor(new.kxs,levels=c(1:nc))

        tab.jxs <- table(new.jxs)
        tab.kxs <- table(new.kxs)

        overlap <- sum(tab.jxs*tab.kxs)
        max.overlap <- length(kxs)*length(jxs)

        cci.null[k,j,i] <- overlap/max.overlap

      }
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
}

