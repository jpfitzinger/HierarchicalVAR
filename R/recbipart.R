# This function recursively bisects the covariance matrix and calculates weights. Arguments include:
#   cov = covariance matrix in raw form (not diagonalised - this is done internally)
#   sortIx = cluster order
#   sortVal = cluster height (this is used to allow for splitting along the hierarchical tree)
#   thresh.val, maxnsplits, const, stop.rule provide additional functionality that is explained later.

recbipart <- function(x, exog, cov, corr, sortIx, sortVal, thresh.val, maxnsplits, stop.rule) {

  # Initialize the output objects
  x.all <- NULL
  sign.all <- NULL
  eq.contains <- list()
  if (thresh.val == 1) thresh.val <- 0.999
  level.counter <- NULL
  level.multiple <- NULL

  # Create recursion function within parent function to avoid use of GlobalEnv
  recbipart.inner <- function(x, level, nlevel, cov, corr, sortIx, sortVal, thresh.val, maxnsplits, stop.rule) {

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
      sign <- rep(0, ncol(x))
      for (j in get(paste0("cItems",i))){
        j0 <- get(paste0("cItems",i))[1]
        if (corr[j0,j] < 0) sign[j] <- -1 else sign[j] <- 1
      }
      sign.all <<- cbind(sign.all, sign)

      if (length(get(paste0("cItems",i)))==1) {
        assign(paste0("var",i), x[,get(paste0("cItems",i))] * (ifelse(level==1,1,grps*level)))
      } else {
        assign(paste0("var",i), rowSums(t(t(x[,get(paste0("cItems",i))] * (ifelse(level==1,1,grps*level))) * sign[get(paste0("cItems",i))])))
      }
      x.grouped[,paste0("var",i)] <- get(paste0("var",i))
    }
    x.all <<- cbind(x.all, as.matrix(x.grouped))
    eq.contains.inner <- c()
    for (i in 1:grps) eq.contains.inner <- c(eq.contains.inner, get(paste0("cItems",i)))
    for (i in 1:grps) eq.contains[[length(eq.contains)+1]] <<- eq.contains.inner

    level.counter <<- c(level.counter, rep(nlevel,grps))
    level.multiple <<- c(level.multiple, rep(ifelse(level==1,1,grps*level), grps))

    for (i in 1:grps) {
      cItems <- get(paste0("cItems",i))
      cItemsVal <- get(paste0("cItemsVal",i))
      if (length(cItems) > stop.rule) {
        recbipart.inner(x, level=grps*level, nlevel=nlevel+1, cov, corr, cItems, cItemsVal, thresh.val, maxnsplits, stop.rule)
      }
      if (length(cItems) <= stop.rule & length(cItems) > 1) {
        recbipart.inner(x, level=grps*level, nlevel=nlevel+1, cov, corr, cItems, cItemsVal, thresh.val = 0, maxnsplits = length(cItems) - 1, stop.rule)
      }
    }

  }

  # run recursion function
  recbipart.inner(x, 1, 1, cov, corr, sortIx, sortVal, thresh.val, maxnsplits, stop.rule)

  x.det <- cbind(x.all, exog)
  XX <- t(x.det) %*% x.det
  XX.temp <- XX[1:ncol(x.all), 1:ncol(x.all)]
  for (i in 1:ncol(x.all)) {
    XX.temp[!sapply(eq.contains, function(x) all(x%in%eq.contains[[i]])),i] <- 0
  }
  XX[1:ncol(x.all), 1:ncol(x.all)] <- XX.temp
  # if (!is.null(exog)) {
  #   if (ncol(exog) > 1) {
  #     XX[-c(1:ncol(x.all)), 1:ncol(x.all)] <- sweep(XX[-c(1:ncol(x.all)), 1:ncol(x.all)], 2, level.multiple, "*")
  #   } else {
  #     XX[-c(1:ncol(x.all)), 1:ncol(x.all)] <- XX[-c(1:ncol(x.all)), 1:ncol(x.all)] * level.multiple
  #   }
  # }

  return(list(XX = XX, sign.all = sign.all, x.all = x.all, level.counter = level.counter, level.multiple = level.multiple))

}
