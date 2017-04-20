#########################################
### functions for single cell RNA seq ###
#########################################
cell_distance <- function(dataframe){
  # calculate cell differences within the defined dataframe
  # assumes all cells will be used, genes in rows, cells in columns
  rho <- cor(dataframe, method="spearman")
  d.rho <- sqrt((1 - rho)/2)
  
  return(d.rho)
}


find_hvg <- function(dataframe, plot=FALSE, p.threshold=1e-2){
  # define a set of highly variable gene for the GFP+ and GFP- separately
  means <- rowMeans(dataframe, na.rm = T)
  vars <- apply(dataframe, 1, var, na.rm=T)
  cv2 <- vars/(means^2)
  
  minMeanForFit <- unname(quantile(means[which(cv2 > 0.2)], 0.8))
  
  # select genes with mean value greater than min value for fitting
  useForFit <- 1/means == 0
  
  # fit with a gamma-distributed GLM
  fit <- glmgam.fit(cbind(a0 = 1, a1tilde=1/means[!useForFit]), cv2[!useForFit])
  
  # calculate % variance explained by the model fit
  resid.var <- var(fitted.values(fit) - cv2[!useForFit])
  total.var <- var(cv2[!useForFit])
  
  # get fitted values and mean-dispersion dependence line
  a0 <- unname(fit$coefficients["a0"])
  a1 <- unname(fit$coefficients["a1tilde"])
  
  xg <- seq(0, max(means[means != Inf]), length.out=100000)
  vfit <- (a1/xg) + a0
  
  # add confidence intervals
  d.f <- ncol(dataframe) - 1
  
  # rank genes by the significance of their deviation from the fit
  # to call HVGs
  a.fit <- (a1/means) + a0
  varFitRatio <- vars/(a.fit * means^2)
  varOrder <- order(varFitRatio, decreasing=T)
  
  oed <- dataframe[varOrder, ]
  
  if(plot == TRUE){
    smoothScatter(means, cv2, xlab="Mean expression", ylab=expression("CV"^2))
    lines(xg, vfit, col="black", lwd=3 )
    lines(xg, vfit * qchisq(0.975, d.f)/d.f, lty=2, col="black")
    lines(xg, vfit * qchisq(0.025, d.f)/d.f,lty=2,col="black")
    # display the 100 most highly variable genes
    points(means[varOrder[1:100]], cv2[varOrder[1:100]], col='red')
  }
  
  pvals <- pchisq(varFitRatio * d.f, d.f, lower.tail = F)
  pvals[is.na(pvals)] <- 1.0
  adj.pvals <- p.adjust(pvals, method='fdr')
  HVG <- adj.pvals <= p.threshold
  return(HVG)
}


create_pseudotime <- function(dataframe, geneset, k, n_eigs, distance, dcs){
  set.seed(24021982)
  dmap <- DiffusionMap(t(dataframe[geneset, ]), n_eigs=n_eigs, k=k, distance=distance)
  dpt <- DPT(dmap)
  dc.map <- as.data.frame(dpt)
  colnames(dc.map) <- paste0("DC", 1:dim(dc.map)[2])
  rownames(dc.map) <- colnames(dataframe)
  
  # reset the seed in case it is used elsewhere
  set.seed(seed=NULL)
  return(dc.map[, 1:dcs])
}


window_cell_distance <- function(dataframe, pseudotime,
                                 direction="reverse"){
  # return a dataframe of windows with the cell-to-cell distances within each window
  n.cells <- dim(dataframe)[2]
  window.size <- round(n.cells * 0.02)
  step.size <- round(window.size/2)
  n.windows <- round(n.cells/step.size)
  
  pseudotime.order <- colnames(dataframe)[order(pseudotime, decreasing=FALSE)]
  
  cell.dist <- list()
  cell.start <- 1
  cell.end <- window.size
  
  for(i in 1:n.windows){
    if((cell.end - step.size) < n.cells){
      window.cells <- pseudotime.order[cell.start:cell.end][!is.na(pseudotime.order[cell.start:cell.end])]
      window.data <- dataframe[, window.cells]
      window.dist <- cell_distance(window.data)
    
      cell.start <- cell.start + step.size
      cell.end <- cell.end + step.size
      cell.dist[[i]] <- window.dist[upper.tri(window.dist, diag=FALSE)]
    }
  }
  
  dist.df <- melt(do.call(cbind, cell.dist))
  colnames(dist.df) <- c("Num", "Window", "Dist")
  return(dist.df)
}


expression_rate <- function(dataframe, pseudotime,
                            tolerance=1e-5,
                            direction="reverse"){
  # return a list of dataframes containing GAM fit, first and second order derivatives
  # pseudotime is a vector of values corresponding to an ordering along a trajectory
  # values should be in the same order as the columns of dataframe
  
  n.cells <- dim(dataframe)[2]
  window.size <- round(n.cells * 0.02)
  step.size <- round(window.size/2)
  n.windows <- round(n.cells/step.size)
  
  pseudotime.order <- colnames(dataframe)[order(pseudotime, decreasing=FALSE)]
  
  n.genes <- dim(dataframe)[1]
  gam.list <- list()
  delta1.list <- list()
  delta2.list <- list()
  genes.av <- list()
  genes.sd <- list()
  cell.dist <- list()
  
  for(k in 1:n.genes){
    gene.data <- dataframe[k, ]
    av.list <- list()
    sd.list <- list()
    cell.start <- 1
    cell.end <- window.size
    
    for(i in 1:n.windows){
      if((cell.end - step.size) < n.cells){
        window.cells <- pseudotime.order[cell.start:cell.end]
        window.data <- as.numeric(gene.data[window.cells])
        window.av <- mean(window.data)
        window.sd <- sd(window.data)
        cell.start <- cell.start + step.size
        cell.end <- cell.end + step.size
        av.list[[i]] <- window.av
        sd.list[[i]] <- window.sd
      }
    }
    
    # fit a GAM to the window
    # print(length(unlist(av.list)))
    # print(1:length(unlist(av.list)))
    gam.data <- data.frame(cbind(c(1:length(unlist(av.list))), unlist(av.list)))
    # gam.data <- cbind(1:n.cells, gam.data[1, ])
    colnames(gam.data) <- c("x", "y")
    window.gam <- gam(y ~ s(x), data=gam.data)
    
    # # code for calculating first and second derivatives from a GAM courtesy of `mnel` on StackOverflow:
    # # http://stackoverflow.com/questions/14207250/determining-derivatives-from-gam-smooth-object
    newDF <- with(gam.data, data.frame(x=x))
    B <- predict(window.gam, newDF, type="response", se.fit=TRUE)
    eps <- 1e-5
    X0 <- predict(window.gam, newDF, type='lpmatrix')
    newDFeps_p <- newDF + eps
    
    X1 <- predict(window.gam, newDFeps_p, type='lpmatrix')
    
    # # finite difference approximation of first derivative
    Xp <- (X1 - X0)/eps
    
    # # first derivative
    fd_d1 <- Xp %*% coef(window.gam)
    
    # # second derivative
    newDFeps_m <- newDF - eps
    X_1 <- predict(window.gam, newDFeps_m, type='lpmatrix')
    Xpp <- (X1 + X_1 - 2*X0)/eps^2
    fd_d2 <- Xpp %*% coef(window.gam)
    
    gam.list[[k]] <- B$fit
    delta1.list[[k]] <- fd_d1[, 1]
    delta2.list[[k]] <- fd_d2[, 1]
    genes.av[[k]] <- unlist(av.list)
    genes.sd[[k]] <- unlist(sd.list)
  }
  
  fits.df <- data.frame(do.call(rbind, gam.list))
  rownames(fits.df) <- rownames(dataframe)
  
  delta1.df <- data.frame(do.call(rbind, delta1.list))
  rownames(delta1.df) <- rownames(dataframe)
  
  delta2.df <- data.frame(do.call(rbind, delta2.list))
  rownames(delta2.df) <- rownames(dataframe)
  
  av.df <- data.frame(do.call(rbind, genes.av))
  rownames(av.df) <- rownames(dataframe)
  
  sd.df <- data.frame(do.call(rbind, genes.sd))
  rownames(sd.df) <- rownames(dataframe)
  

  out.list <- list("fit"=fits.df,
                   "derive1"=delta1.df,
                   "derive2"=delta2.df,
                   "average"=av.df,
                   "sigma"=sd.df)
  return(out.list)
}


.temp_cor <- function(x, y){
  n.time <- length(x)
  sum_prod <- c()
  sum_xsq <- c()
  sum_ysq <- c()
  
  xi.v <- diff(as.numeric(x))
  yi.v <- diff(as.numeric(y))
  sum_prod <- sum(xi.v * yi.v)
  sum_xsq <- sum(xi.v ** 2)
  sum_ysq <- sum(yi.v ** 2)
  
  nume <- sum(sum_prod)
  denom <- sqrt(sum(sum_xsq)) * sqrt(sum(sum_ysq))
  
  return(nume/denom)
}

.apply_tempcor <- function(X, dataframe){
  atcor <- apply(dataframe, 1,
                 FUN=function(Q) .temp_cor(Q, X))
  return(atcor)
}

temp_cor <- function(dataframe){
  tcor <- apply(dataframe, 1,
                FUN=function(X) .apply_tempcor(X, dataframe))
  return(tcor)
}
