get.clusters <- function(corr, clust.method = "AGNES", clust.control, absolute = F) {

  if (absolute) corr <- abs(corr)
  distmat <- ((1 - corr) / 2)^0.5
  if (clust.method == "AGNES") {
    if (is.null(clust.control)) {
      clust.control <- list(method = "single")
    }
    clust <- do.call(cluster::agnes, c(list(x = dist(distmat)), clust.control))
    #clust$height <- recTree(corr, clust$order)
  }
  if (clust.method == "DIANA") {
    clust <- do.call(cluster::diana, c(list(x = dist(distmat)), clust.control))
    #clust$height <- recTree(corr, clust$order)
  }
  if (clust.method == "OPT") {
    W <- function(x) {
      corr_o <- corr[x,x]
      dist_o <- ((1-corr_o)/2)^0.5
      i <- matrix(c(1:ncol(corr)),ncol(corr),ncol(corr))
      j <- t(i)
      weight <- (dist_o) * ((i - j)^2) # should be abs corr

      sum(weight)
    }
    tmp <- gaoptim::GAPerm(W,ncol(corr))
    if (is.null(clust.control)) clust.control <- list(niter = 500)
    tmp$evolve(clust.control$niter)
    order <- tmp$bestIndividual()
    clust <- list(order = order, height = recTree(corr,order))
  }

  order <- clust$order
  height <- clust$height
  height <- c(height, height[length(height)])
  list(order = order, height = height)

}

recTree <- function(corr, order) {

  h <- rep(ncol(corr)*10, ncol(corr))

  recTreeInner <- function(sub_idx) {

    corr_o <- corr[sub_idx,sub_idx]
    mac <- sapply(1:(length(sub_idx)-1), function(x) mean(abs(corr_o[1:x,-c(1:x)])))
    split <- c(mac == min(mac),F)
    h[-sub_idx[split]] <<-  h[-sub_idx[split]] - length(sub_idx)

    new_idx <- 1:c(1:length(sub_idx))[split]
    sub_idx1 <- sub_idx[new_idx]
    sub_idx2 <- sub_idx[-new_idx]

    if (length(sub_idx1)>1) recTreeInner(sub_idx1)
    if (length(sub_idx2)>1) recTreeInner(sub_idx2)

  }

  recTreeInner(order)
  h <- h[order]
  return(h[-length(h)])

}



