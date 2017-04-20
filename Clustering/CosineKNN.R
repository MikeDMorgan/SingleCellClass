###################################
## crisp and fuzzy kNN functions ##
###################################

uix <- function(x, xj, m){
  # shouldn't this be the same distance metric used for
  # the original kNN?
  dist <- 1 - cosine(x, xj)
  pow <- (m - 1)
  
  # what is the consequence for varying values of m?
  ui <- 1/(dist ** pow)
  # ui <- 1/dist
  return(ui)
}


fuzzy_cos <- function(train, testInstance, testNN, m){

  class.train <- train[, unique(testNN)]
  if(is.null(dim(class.train))){
    # if there is a single class then this is the
    # equivalent of a crisp kNN or membership = 1.0
    uis <- as.numeric(uix(testInstance, class.train, m))
    names(uis) <- unique(testNN)
  }
  else{
    uis <- apply(class.train, 2,
                 FUN=function(Q) uix(testInstance, Q, m))
  }
  ui_sum <- uis/sum(uis)
  return(ui_sum)
}


assign_classes <- function(test, y.nn, class.average, m){
  # for each cell calculate it's fuzzy membership to
  # the class average
  col.indx <- seq_len(dim(test)[2])
  class_fuzz <- lapply(col.indx,
                       function(Q) fuzzy_cos(class.average, test[, Q],
                                             y.nn[, Q], m))
  return(class_fuzz)  
  
}


class_average <- function(train, cl){
  names(cl) <- colnames(train)
  av.list <- list()
  labels <- unique(cl)
  for(l in 1:length(labels)){
    lbl <- labels[l]
    xprs <- train[, which(cl == lbl, arr.ind=T)]
    if(is.null(dim(xprs))){
      proto <- xprs
    }
    else{
      proto <- rowMeans(xprs)
    }
    
    av.list[[lbl]] <- proto
  }
  av.mat <- as.data.frame(do.call(cbind, av.list))
  
  colnames(av.mat) <- labels
  rownames(av.mat) <- rownames(train)
  return(av.mat)
}


fuzzyCosKNN <- function(train, test, cl, k, m){
  # calculate class averages
  class.av <- class_average(train, cl)
  # enforce a matrix on nn.matrix?
  # get the labels for nearest neighbours
  nn.matrix <- apply(test, 2,
                     FUN=function(P) cosine_knn(train, P, k))
  nn.labels <- matrix(cl[nn.matrix], nrow=k)
  colnames(nn.labels) <- colnames(test)
  
  # calculate fuzzy set membership
  fuzzy_mem <- assign_classes(test, nn.labels, class.av, m)
  
  # pad out non-membership classes as zeros
  all.present <- unique(as.vector(nn.labels))
  n <- length(all.present)
  p <- dim(test)[2]
  empty <- data.frame(matrix(numeric(n*p),
                             ncol=p, nrow=n))
  rownames(empty) <- all.present
  colnames(empty) <- colnames(test)
  
  for(c in seq_len(length(fuzzy_mem))){
    c.k <- fuzzy_mem[c]
    lbls <- names(c.k[[1]])
    empty[lbls, c] <- c.k[[1]]
  }
  return(empty)
  
}


cosine_knn <- function(train, testInstance, k){
  require(lsa)
  # this is not fast
  distance <- apply(train, 2, FUN=function(D) cosine(D, testInstance))

  sort_dist <- sort(1 - distance, decreasing=F)
  nNames <- names(sort_dist)[1:k]
  # need to return indices in train
  # to retrieve their proper class label
  nn <- which(colnames(train) %in% nNames, arr.ind=T)
  
  return(nn)
}


vote <- function(nn.index, labels){
  # nn.index is the index for the k-NN
  # return the label with the most votes
  res <- table(labels[nn.index])
  max.res <- max(res)
  decision <- names(res)[res == max.res]
  
  if(length(decision) == 1){
    return(decision)
  }
  else {
    # randomly select one?
    ran <- floor(runif(1, min=1, max=length(decision)))
    return(decision[ran])
  }
}


cosineKNN <- function(train, test, cl, k){
  # enforce a matrix on nn.matrix?
  nn.matrix <- apply(test, 2,
                     FUN=function(P) cosine_knn(train, P, k))
  kNN <- apply(nn.matrix, 2,
               function(V) vote(V, cl))
  return(kNN)
}


##############################
## Distance-based functions ##
##############################

euclidDist <- function(x, y){
  sqrt(sum((x - y) ** 2))
}

cosine_distance <- function(dataframe){
  # calculate all pairwise cosine distances on a matrix/dataframe
  # calculate column-column
  cos.list <- list()
  for(x in 1:dim(dataframe)[2]){
    col.vec <- dataframe[, x]
    print(length(col.vec))
    cos.dist <- apply(data.frame, 2,
                      FUN=function(D) cosine(D, col.vec))
    cos.list[[x]] <- cos.dist
  }
  out.mat <- as.data.frame(do.call(rbind, cos.list))
  colnames(out.mat) <- colnames(dataframe)
  rownames(out.mat) <- rownames(dataframe)
  return(out.mat)
}


proportional_distance <- function(train, testInstance,
                                  proportional=TRUE){
  # calculate the proportional distance between each
  # class and each test case
  
  distance <- apply(train, 2,
                    FUN=function(D) cosine(D, testInstance))
  if(proportional){
    prop.dist <- distance/sum(distance)
  }
  else{
    prop.dist <- distance
  }
  return(prop.dist)
}


proportionalDistance <- function(train, test, cl,
                                 return.prop=TRUE){
  prop.mat <- apply(test, 2,
                      FUN=function(P) proportional_distance(train, P,
                                                            return.prop))

  rownames(prop.mat) <- cl
  return(prop.mat)
}


####################################
## Equality/specificity functions ##
####################################
Jaccard <- function(dataframe, binary=TRUE){
  # apply the Jaccard coefficient to each pair of
  # rows in a matrix
  # sanity check that diagonals == 1
  
  # there is a linear algebraic calculation for Jaccard
  # XTX = A
  # (notX)T(notX) = D
  # Jaccard = A/(n-D)
  # two transposes, then some element-wise operations.
  if(binary == TRUE){
    dat <- as.matrix(dataframe)
    n <- dim(dataframe)[1]
    A <- t(dat) %*% dat
    notX <- abs(dat - 1)
    D <- t(notX) %*% notX
    
    jac.df <- data.frame(A/(n-D))
  }
  else{
    n.rows <- dim(dataframe)[1]
    jac.df <- data.frame(apply(dataframe, 1,
                               FUN=function(Q) row_jaccard(x=dataframe,
                                                           row=Q,
                                                           binary=binary)))
  }
  return(jac.df)
}


row_jaccard <- function(x, row, binary=TRUE){
  jac.row <- apply(x, 1,
                   FUN=function(C) jaccard_coef(x=C,
                                                y=row,
                                                binary=binary))
  return(jac.row)
}

jaccard_coef <- function(x, y, binary=TRUE){
  # check if the data are binary vectors
  # if not then assume category labels so use set operations
  
  if(binary == TRUE){
    numerator <- sum(x == 1 & y == 1)
    denom <- sum(sum(x == 1 & y == 1),
                 sum(x == 1 & y == 0),
                 sum(x == 0 & y == 1))
  }
  else{
    numerator <- length(intersect(x, y))
    denom <- sum(length(intersect(x, y)),
                 length(setdiff(x, y)),
                 length(setdiff(y, x)))
    
  }
  jac <- numerator/denom
  return(jac)
}

Gini <- function(x){
  # absolute differences between each pair of points
  diff.vec <- numeric(0)
  for(i in 1:length(x)){
    for(j in 1:length(x)){
      if(i != j){
        diff.vec <- c(diff.vec, abs(x[i] - x[j]))
      }
    }
  }
  sum.x <- sum(x)
  n <- length(x)
  gini.c <- sum(diff.vec)/((2*n) * sum.x)
  return(gini.c)
}


Tau <- function(x){
  # calculate specificity index over vector/expression profile
  max_x <- max(x)
  max_I <- 1 - (x/max_x)
  numerator <- sum(max_I)
  denominator <- length(x)
  
  tau <- numerator/denominator
  
  return(tau)
}


#############################
## Normalisation functions ##
#############################

euclid_norm <- function(vector){
  # calculate the Euclidean norm scaled vector
  norm <- sqrt(t(vector) %*% vector)
  return(vector/norm)
}


binary_convert <- function(x, threshold){
  return(as.integer(x >= threshold))
}

ternary_convert <- function(x, threshold){
  bx <- x
  bx[x <= threshold] <- -1
  bx[is.na(x)] <- 0
  bx[x > threshold] <- 1
  return(bx)
}

#################################
## Gene contribution functions ##
#################################

geneContributions <- function(input.matrix, train, test, cl){
  # matrix is classes X cells
  # non-zero entries represent fuzzy assignment
  # for that cell
  
  # train should be the class average, or calculate them
  # on the fly?
  class.av <- class_average(train, cl)
  classes <- rownames(input.matrix)
  cell.vec <- seq_len(dim(test)[2])
  # turn input.matrix into binary matrix
  bin.mat <- matrix(as.integer(input.matrix != 0),
                    ncol=dim(input.matrix)[2],
                    nrow=dim(input.matrix)[1])
  
  sapply(cell.vec, 
         FUN=function(H) class_genes(bin.mat[, H],
                                     test[, H],
                                     class.av,
                                     classes))
  
}


class_genes <- function(bin.vec, testInstance, train, cl){
  test.unit <- euclid_norm(testInstance)
  classes <- cl[as.logical(bin.vec)]
  class.list <- list()
  for(cls in classes){
    train.unit <- euclid_norm(train[, cls])
    
    # either find a better way to set this threshold,
    # or allow the user to set it
    test.thresh <- as.integer(test.unit > 0)
    train.thresh <- as.integer(train.unit > 0)
    
    unit.mat <- test.thresh %*% t(train.thresh)
    mat.diag <- diag(unit.mat)
    exp.genes <- unique(rownames(train)[as.logical(mat.diag)])
    class.list[[cls]] <- exp.genes
  }
  return(class.list)
}

grepREACT <- function(id, mapkeys){
  unique(unlist(mapkeys[id], use.names=F))
}


mapGene2Reactome <- function(...){
  require(goseq)
  require(reactome.db)
  require(org.Mm.eg.db)
  
  # map ensembl:entrez, then entrez:reactome
  # need to generalise this to other species
  en2eg <- as.list(org.Mm.egENSEMBL2EG)
  eg2reactome <- as.list(reactomeEXTID2PATHID)
  react <- lapply(en2eg, grepREACT, eg2reactome)
  return(react)
  
}


cellClassGO <- function(gene.list, back.genes,
                        gene.lengths){
  # reverse entrez:reactome mapping
  # only do once, it takes ages otherwise
  react <- mapGene2Reactome()
  cell.list <- list()
  for(cell in seq_len(length(gene.list))){
    enrich.list <- lapply(gene.list[[cell]],
                          FUN=function(T) geneClassGO(T,
                                                      back.genes,
                                                      gene.lengths,
                                                      react))
    cell.list[[cell]] <- enrich.list
  }
  return(cell.list)
}


geneClassGO <- function(test.genes, back.genes,
                        gene.lengths, mapping){
  # required extensions:
  # different species, different gene ID mappings,
  # retrieve gene lengths from txDB object
  # perform enrichment testing of reactome terms with goseq
  gene.vector <- as.integer(unique(back.genes) %in% unique(test.genes))
  names(gene.vector) <- unique(back.genes)
  all.gene.vector <- gene.vector[names(gene.lengths)]
  
  pwf.all <- nullp(all.gene.vector, "mm10",
                   "ensGene", bias.data=gene.length,
                   plot.fit=FALSE)
  GO.wall.all <- suppressMessages(goseq(pwf.all, genome="mm10",
                                        id="ensGene", gene2cat=mapping,
                                        method="Wallenius", use_genes_without_cat = T))
  
  GO.wall.all$padjust <- p.adjust(GO.wall.all$over_represented_pvalue, method="BH")
  GO.wall.all$foldEnrich <- ((GO.wall.all$numDEInCat/length(test.genes))/(GO.wall.all$numInCat/length(gene.length)))
  GO.sig.all <- GO.wall.all[(GO.wall.all$padjust <= 0.05), ]
  
  id2name <- as.list(reactomePATHID2NAME)
  GO.sig.all$description <- do.call(rbind, id2name[GO.sig.all$category])
  
  return(GO.sig.all)
}

####################################################
## TF enrichment testing in tissue enriched mTECs ##
####################################################
test_geneset_overlap <- function(geneset, genes.list, tissue.genes,
                                 tissue.bin){
  geneset.store <- list()
  for(i in 1:length(genes.list)){
    if(length(genes.list[[i]]) > 1){
      cell.tis <- names(genes.list[[i]])
      for(tis in 1:length(cell.tis)){
        en.tis <- cell.tis[[tis]]
        shared.genes <- genes.list[[i]][[en.tis]]
        tis.genes <- tissue.genes[as.logical(tra.bin[, en.tis])]
        cell.geneset <- intersect(shared.genes, geneset)
        tis.geneset <- intersect(tis.genes, geneset)
        geneset.store[["Multiple"]] <- c(geneset.store[["Multiple"]],
                                        length(cell.geneset)/length(tis.geneset))
      }
    }
  else{
    cell.tis <- names(genes.list[[i]])
    shared.genes <- genes.list[[i]][[cell.tis]]
    tis.genes <- tissue.genes[as.logical(tra.bin[, cell.tis])]
    cell.geneset <- intersect(shared.genes, geneset)
    tis.geneset <- intersect(tis.genes, geneset)
    geneset.store[["Single"]] <- c(geneset.store[["Single"]],
                                  length(cell.geneset)/length(tis.geneset))
  }
  }
  return(geneset.store)
}


tf.fisher_test <- function(tf.genes, tissue.genes, cell.genes, all.genes){
  # test whether there is an enrichment/depletion of TF genes
  # in the tissue-specific genes in each mTEC:tissue pair
  Mg <- cell.genes
  Tg <- tf.genes
  Pg <- tissue.genes
  
  # calculate 2x2 contigency table
  aF <- length(intersect(cell.genes, tf.genes))
  bF <- length(setdiff(cell.genes, tf.genes))
  cF <- length(setdiff(tf.genes, cell.genes))
  dF <- length(setdiff(setdiff(all.genes, cell.genes), tf.genes))
  fish <- matrix(c(aF, bF, cF, dF), nrow=2, ncol=2)
  fish.res <- fisher.test(fish)
  fish.p <- fish.res$p.value
  fish.or <- fish.res$estimate
  fish.ci <- fish.res$conf.int
  fish.vec <- c(fish.or, fish.ci, fish.p)
  names(fish.vec) <- c("OR", "LCI", "UCI", "P")
  return(fish.vec)
}

adjacency.matrix <- function(geneset.1, geneset.2){
  # calculate an adjacency matrix as the proportion of times
  # two genes are found together.  If there is only one group of
  # genes then set geneset.1 and geneset.2 to the same matrix
  # of genes X cells
  n.cells <- dim(geneset.1)[2]
  gene.list <- list()
  for(k in 1:(dim(geneset.1)[1])){
    vec <- geneset.1[k, ]
    gene.list[[k]] <- rowSums(geneset.2 * t(vec))
  }
  gene.mat <- as.data.frame(do.call(rbind, gene.list))/n.cells
  rownames(gene.mat) <- rownames(geneset.1)
  return(gene.mat)
}