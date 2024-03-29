
#' @export
#'
#' @title Default algorithm for cluster cohesion index
#'
#' @param V The scGBM scores
#' @param nstart A \code{kmeans} parameter
#' @param centers A \code{kmeans} parameter
#'
#' @return The function applies \code{kmeans} to \code{V} and returns
#' the assignments of clusters
#'
#' @author Phillip B. Nicol <philnicol740@gmail.com>
cluster.default <- function(V,nstart=25,centers=5) {
  kmeans(V,nstart=nstart,centers=centers)$cluster
}

#' @export
#'
#' @title Evaluate cluster stability using the Cluster Cohesion Index (CCI)
#'
#' @param gbm The scGBM object estimated from \code{gbm.sc}.
#' @param cluster.orig A vector of length J containing the labeling of
#' cells into clusters.
#' @param reps The number of simulated datasets to make.
#' @param cluster.fn A function that takes as input a set of (simulated) scores V and
#' returns a labeling of cells into clusters.
#'
#' @return A list with components
#' \itemize{
#' \item \code{H.table} - The values used to make the cluster cohesion
#' heatmap.
#' \item \code{cci_diagonal} - A ggplot object that plots the distribution of intra-cluster
#' CCI.
#' \item \code{cci} - A 3-dimensional array of all CCI values (inter-and-intra-cluster).
#' \item \code{coarse_cluster} - A vector of length J containing the new coarse cluster
#' labels where clusters with high inter-CCI have been combined into a single cluster.
#' Disclaimer: this part is experimental and has not been thoroughly benchmarked.
#' }
#'
#' @author Phillip B. Nicol <philnicol740@gmail.com>
CCI <- function(gbm,
                cluster.orig,
                reps=100,
                cluster.fn=cluster.default,
                null.dist = FALSE,
                ...) {
  nc <- length(table(cluster.orig))
  cci <- array(0,dim=c(nc,nc,reps))
  H.table <- matrix(0,nrow=nc,ncol=nc)
  V <- gbm$scores
  se_V <- gbm$se_scores
  J <- nrow(V)
  M <- ncol(V)

  if(null.dist) {
    V.null <- matrix(rnorm(n=J*M,sd=se_V),nrow=J,ncol=M)
    rownames(V.null) <- rownames(V); colnames(V.null) <- colnames(V)

    cluster.null <- cluster.fn(V.null,...)

    nc.null <- length(table(cluster.null))
    cci.null <- array(0,dim=c(nc.null,nc.null,reps))
  }

  for(i in 1:reps) {
    cat("Iteration ", i, "\n")
    M <- length(gbm$D)
    e <- matrix(rnorm(n=J*M,sd=se_V),nrow=J,ncol=M)

    V.new <- V+e
    cluster <- cluster.fn(V.new,...)
    nc.new <- length(unique(cluster)) #New number of clusters

    for(k in 1:nc) {
      for(j in 1:nc) {
        kxs <- which(cluster.orig==unique(cluster.orig)[k])
        jxs <- which(cluster.orig==unique(cluster.orig)[j])

        new.jxs <- cluster[jxs]
        new.kxs <- cluster[kxs]

        new.jxs <- factor(new.jxs,levels=c(1:nc.new))
        new.kxs <- factor(new.kxs,levels=c(1:nc.new))

        tab.jxs <- table(new.jxs)
        tab.kxs <- table(new.kxs)

        overlap <- sum(tab.jxs*tab.kxs)
        max.overlap <- length(kxs)*length(jxs)

        cci[k,j,i] <- overlap/max.overlap

      }
    }

    if(null.dist) {
      e <- matrix(rnorm(n=J*M,sd=se_V),nrow=J,ncol=M)
      V.new <- V.null + e
      cluster <- cluster.fn(V.new,...)
      nc.new <- length(unique(cluster)) #New number of clusters

      for(k in 1:nc.null) {
        for(j in 1:nc.null) {
          kxs <- which(cluster.null==unique(cluster.null)[k])
          jxs <- which(cluster.null==unique(cluster.null)[j])

          new.jxs <- cluster[jxs]
          new.kxs <- cluster[kxs]

          new.jxs <- factor(new.jxs,levels=c(1:nc.new))
          new.kxs <- factor(new.kxs,levels=c(1:nc.new))

          tab.jxs <- table(new.jxs)
          tab.kxs <- table(new.kxs)

          overlap <- sum(tab.jxs*tab.kxs)
          max.overlap <- length(kxs)*length(jxs)

          cci.null[k,j,i] <- overlap/max.overlap
        }
      }
    }
  }

  #Create heatmap
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
  p <- p + labs(fill="Mean")

  rownames(H.table) <- unique(cluster.orig)
  colnames(H.table) <- unique(cluster.orig)
  pheatmap::pheatmap(H.table,legend=TRUE, color=colorRampPalette(c("white","red"))(100),
           breaks=seq(0,1,by=0.01),
           rownames=TRUE,
           colnames=TRUE,
           cluster_rows = ifelse(nrow(H.table) > 1, TRUE,FALSE),
           cluster_cols = ifelse(nrow(H.table) > 1, TRUE,FALSE))

  #Standard CCI
  cci.orig <- matrix(0,nrow=reps,ncol=nc)
  for(i in 1:nc) {
    cci.orig[,i] <- cci[i,i,]
  }

  if(null.dist) {
    cci.null95 <- c()
    for(i in 1:nc.null) {
      cci.null95 <- c(cci.null95, cci.null[i,i,])
    }

    cci.null95 <- quantile(cci.null95, 0.95)
  }

  out <- list()
  #out$heatmap <- p
  out$H.table <- H.table

  coarse_cluster <- rep(0,J)
  ctr <- 1
  H.cut <- as.matrix(H.table)
  #print(H.cut)
  diag(H.cut) <- 0
  H.cut[H.cut < 0.5/nc] <- 0
  G <- igraph::graph_from_adjacency_matrix(H.cut,
                                           mode="undirected",
                                           weighted=TRUE)
  while(length(igraph::V(G)) > 0) {
    vs <- igraph::largest.cliques(G)[[1]]
    coarse_cluster[cluster.orig %in% names(vs)] <- ctr
    #print(table(coarse_cluster[cluster.orig %in% vs]))
    ctr <- ctr + 1
    G <- igraph::delete.vertices(G,vs)
  }

  ##Plotting
  colnames(cci.orig) <- unique(cluster.orig)
  df <- reshape2::melt(cci.orig)
  df$Var2 <- as.factor(df$Var2)

  p <- ggplot(data=df,aes(x=Var2,y=value,fill=Var2))
  p <- p + geom_boxplot(alpha=0.5)
  p <- p + theme_bw()
  if(null.dist) {
    p <- p + geom_hline(yintercept=cci.null95, linetype="dashed",color="blue")
  }
  p <- p + xlab("Cluster") + ylab("CCI")
  p <- p + ggtitle("Distribution of Cluster Cohesion Index")
  p  <- p + guides(fill="none")
  p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  out$cci_diagonal <- p
  out$cci <- cci
  if(null.dist) {
    out$cci.null <- cci.null
    out$cci.null95 <- cci.null95
    out$cluster.null <- cluster.null
  }
  out$coarse_cluster <- coarse_cluster
  return(out)
}
