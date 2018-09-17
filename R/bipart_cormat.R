bipart.cormat <- function(sortIx, sortVal, maxnsplits, thresh.val, min.cluster) {

  # Convert "threshold" from 0-1 value to units of cluster height
  interval <- sortVal[-length(sortVal)]
  thresh <- min(interval) + thresh.val*(max(interval) - min(interval))

  # Find qualifying splits (based on threshold)
  grps <- sum(interval>=thresh)+1
  splits <- c(1:length(sortIx))[c(interval>=thresh,F)]

  # Find splits that maximise diversification
  #   If there are more qualifying splits than maxnsplits, the goal is to pick those that maximise tree symmetry
  if (maxnsplits < length(splits)) {
    split.combo <- combn(splits, maxnsplits, FUN = sort)
    split.augment <- rbind(0, split.combo, length(sortIx))
    split.size <- split.augment[-1,] - split.augment[-nrow(split.augment),]
    split.diverse <- colMeans(split.size^2)
    if (maxnsplits == 1) {
      opt.splits <- split.combo[split.diverse==min(split.diverse)]
    } else {
      opt.splits <- split.combo[,split.diverse==min(split.diverse)] %>% as.vector
    }
    grps <- maxnsplits + 1
    splits <- opt.splits
  }

  # Generate sort and value vectors for each sub-group
  #   This uses the function "getClusterVar" which is described below
  x.grouped <- data.frame(matrix(NA, nrow(x), grps))
  colnames(x.grouped) <- paste0("var",1:grps)
  for (i in 1:grps) {
    st <- if (i==1) {1} else {splits[i-1]+1}
    en <- if (i==grps) {length(sortIx)} else {splits[i]}
    range <- st:en
    assign(paste0("subIdx",i), range)
    assign(paste0("cItems",i), sortIx[get(paste0("subIdx",i))])
    assign(paste0("cItemsVal",i), sortVal[get(paste0("subIdx",i))])

    # Sign adjustment
    for (j in get(paste0("cItems",i))){
      j0 <- get(paste0("cItems",i))[1]
      if (corr[j0,j] < 0) sign[j] <<- -1 else sign[j] <- 1
    }

    if (length(get(paste0("cItems",i)))==1) {
      assign(paste0("var",i), x[,get(paste0("cItems",i))])
    } else {
      assign(paste0("var",i), rowSums(t(t(x[,get(paste0("cItems",i))]) * sign[get(paste0("cItems",i))])))
    }
    x.grouped[,paste0("var",i)] <- get(paste0("var",i))
  }


}
