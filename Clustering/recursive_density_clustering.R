iterDensityCluster <- function(dist.matrix, delta, rho){
  dens <- densityClust(dist.matrix)
  rho_t <- rho * 0.01
  delta_t <- delta * 0.25
  dens <- findClusters(dens,
                       delta=delta_t,
                       rho_t)
  clusters <- dens$clusters
  
  return(dens)
}

recursiveDensityClustering <- function(input.matrix, limit=10){
  # recursively cluster data points, dropping the largest cluster
  # on each iteration
  
  # setup the initial clustering
  # how to select the initial values of delta and rho?
  # 5% of the total highest value?
  # lots of clusters to start to maximise granularity
  
  require(densityClust)
  dist.matrix <- as.dist(1 - input.matrix)
  # set initial delta values
  init.dens <- densityClust(dist.matrix)
  init.delta <- max(init.dens$delta)
  init.rho <- max(init.dens$rho)
  
  dens <- iterDensityCluster(dist.matrix,
                             delta=init.delta,
                             rho=init.rho)
  gamma.v <- dens$delta * dens$rho
  
  # cluster centers have greatest density and distance, i.e. gamma > inf
  log.gamma <- na.omit(log(gamma.v[order(gamma.v, decreasing=T)]))
  log.nz.gamma <- log.gamma[log.gamma > -1000]
  
  # big cluster centres fall below the mean log(gamma) + 2xsigma
  # order gamma values, preserve names of ordering
  # cluster centers should lie out on the tails of the distribution
  c.dens <- names(log.nz.gamma[log.nz.gamma >= mean(log.nz.gamma) + 2*sqrt(var(log.nz.gamma))])
  cluster_match <- as.data.frame(cbind(colnames(input.matrix), dens$clusters))
  rownames(cluster_match) <- cluster_match$V1
  top.clust <- as.character(unique(cluster_match[intersect(c.dens, cluster_match$V1), ]$V2))
  print(top.clust)
  cluster.list <- list()
  cluster_match$clusters <- sapply(cluster_match$V2,
                                   FUN=function(X) paste(1, X, sep='.'))
  cluster.list[[1]] <- cluster_match[cluster_match$V2 == top.clust, ]
  
  new.input <- input.matrix[which(cluster_match$V2 != top.clust), 
                            which(cluster_match$V2 != top.clust)]
  rownames(new.input) <- colnames(new.input)
  new.dist <- as.dist(1 - new.input)
  plot(dens$rho, dens$delta)
  abline(v=rep(max(dens$rho)*0.01, 80), lty=3, lwd=2, col='purple')
  lines(rep(max(dens$delta)*0.1, 80), lty=3, col='red', lwd=2)
  
  
  for(k in 2:limit){
    # drop all clusters with value above gamma_threshold?
    new.delta <- max(dens$delta)
    new.rho <- max(dens$rho)
    dens <- iterDensityCluster(new.dist,
                               delta=new.delta,
                               rho=new.rho)
    
    gamma.v <- dens$delta * dens$rho
    
    # cluster centers have greatest density and distance, i.e. gamma > inf
    log.gamma <- na.omit(log(gamma.v[order(gamma.v, decreasing=T)]))
    log.nz.gamma <- log.gamma[log.gamma > -1000]
    
    # big cluster centres fall below the mean log(gamma) + 2xsigma
    # order gamma values, preserve names of ordering
    # cluster centers should lie out on the tails of the distribution
    c.dens <- names(log.nz.gamma[log.nz.gamma >= mean(log.nz.gamma) + 2*sqrt(var(log.nz.gamma))])
    cluster_match <- as.data.frame(cbind(colnames(new.input), dens$clusters))
    rownames(cluster_match) <- cluster_match$V1
    top.clust <- as.character(unique(cluster_match[intersect(c.dens, cluster_match$V1), ]$V2))
    print(top.clust)
    cluster_match$clusters <- sapply(cluster_match$V2,
                                    FUN=function(X) paste(k, X, sep='.'))
    cluster.list[[k]] <- cluster_match[cluster_match$V2 == top.clust, ]
    
    plot(dens$rho, dens$delta)
    abline(v=rep(max(dens$rho)*0.01, 80), lty=3, lwd=2, col='purple')
    lines(rep(max(dens$delta)*0.1, 80), lty=3, col='red', lwd=2)
    
    new.mat <- new.input[which(cluster_match$V2 != top.clust),
                         which(cluster_match$V2 != top.clust)]
    colnames(new.mat) <- rownames(cluster_match)[which(cluster_match$V2 != top.clust)]
    rownames(new.mat) <- rownames(cluster_match)[which(cluster_match$V2 != top.clust)]
    print(dim(new.mat))
    new.input <- new.mat
    new.dist <- as.dist(1 - new.input)
    
  }
  cluster.list[[k+1]] <- cluster_match[cluster_match$V2 != top.clust, ]
  return(cluster.list)
}

par(mar=c(4, 3, 2, 2))
par(mfrow=c(3, 2))
limits <- 20
test_dist <- var1000.cor
clusters_list <- recursiveDensityClustering(test_dist, limit=limits)
iter.clust <- do.call(rbind, clusters_list)
print(randIndex(table(iter.clust$clusters,
                      level1_df$variable),
                adjust=T))
# pheatmap(table(iter.clust$clusters,
#                level1_df$variable)/rowSums(table(iter.clust$clusters,
#                                                                      level1_df$variable)))