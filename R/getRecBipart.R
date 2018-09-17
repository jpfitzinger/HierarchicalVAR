# This function recursively bisects the covariance matrix and calculates weights. Arguments include:
#   cov = covariance matrix in raw form (not diagonalised - this is done internally)
#   sortIx = cluster order
#   sortVal = cluster height (this is used to allow for splitting along the hierarchical tree)
#   thresh.val, maxnsplits, const, stop.rule provide additional functionality that is explained later.

getRecBipart <- function(y, x, exog, cov, corr, sortIx, sortVal, thresh.val, maxnsplits, stop.rule, p.crit = 0.05, absolute = T, det = "before") {

  # Initialize the coefficient vector to 0
  dim.x <- ncol(x) + ifelse(!is.null(exog), ncol(exog), 0)
  w <- rep(0, dim.x)
  sign <- rep(1,ncol(x))
  if (thresh.val == 1) thresh.val <- 0.999

  if (!is.null(exog) & det == "before") {
    mod.exog <- lm(y~ -1 + exog)
    y <- residuals(mod.exog)
    w[-c(1:ncol(x))] <- coefficients(mod.exog)
  }

  # Create recursion function within parent function to avoid use of GlobalEnv
  recurFun <- function(y, x, exog, cov, corr, sortIx, sortVal, thresh.val, maxnsplits, stop.rule, p.crit, det) {

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
        if (corr[j0,j] < 0) sign[j] <<- -1 else sign[j] <<- 1
      }

      if (length(get(paste0("cItems",i)))==1) {
        assign(paste0("var",i), x[,get(paste0("cItems",i))])
      } else {
        assign(paste0("var",i), rowSums(t(t(x[,get(paste0("cItems",i))]) * sign[get(paste0("cItems",i))])))
      }
      x.grouped[,paste0("var",i)] <- get(paste0("var",i))
    }

    # Calculate alpha
    if (!is.null(exog) & det=="during") data <- as.data.frame(cbind(y, x.grouped, exog)) else data <- as.data.frame(cbind(y, x.grouped))
    mod <- lm(as.formula("y~ -1 + ."), data = data)
    coefs <- coefficients(mod)

    pvals <- summary(mod)$coefficients[,4]
    pvals <- pvals <= p.crit
    coefs[!pvals] <- 0

    # Apply weights - assign up a level into the parent environment
    # if (!is.null(exog)) w[-c(1:ncol(x))] <<- w[-c(1:ncol(x))] + coefs[-c(1:ncol(x.grouped))]
    coefs.exog <- coefs[-c(1:ncol(x.grouped))]
    coefs <- coefs[1:ncol(x.grouped)]
    w.temp <- w[1:ncol(x)]
    for (i in 1:grps) {
      w.temp[get(paste0("cItems",i))] <- w.temp[get(paste0("cItems",i))] + rep(coefs[i], length(get(paste0("cItems",i)))) * sign[get(paste0("cItems",i))]
    }
    w[1:ncol(x)] <<- w.temp

    if (!is.null(exog) & det=="during") {
      y <- residuals(mod) + exog %*% coefs.exog
      #y <- residuals(mod)
    } else {
      y <- residuals(mod)
    }
    #y <- y / grps

    # Calculate group internal variation
    w.a <- vector(mode = "logical", length = grps)
    for (i in 1:grps) {
      cItems <- get(paste0("cItems",i))
      corr.a <- corr[cItems,cItems]
      w.a[i] <- mean(abs(corr.a))
      w.a[i] <- 1
    }
    w.a <- (1/w.a) / sum(1/w.a)

    for (i in 1:grps) {
      cItems <- get(paste0("cItems",i))
      cItemsVal <- get(paste0("cItemsVal",i))
      if (length(cItems) > stop.rule) {
        recurFun(y*w.a[i], x, exog = exog, cov, corr, cItems, cItemsVal, thresh.val, maxnsplits, stop.rule, p.crit, det)
      }
      if (length(cItems) <= stop.rule & length(cItems) > 1) {
        recurFun(y*w.a[i], x, exog = exog, cov, corr, cItems, cItemsVal, thresh.val = 0, maxnsplits = length(cItems) - 1, stop.rule, p.crit, det)
      }
    }

  }

  # run recursion function
  recurFun(y, x, exog = exog, cov, corr, sortIx, sortVal, thresh.val, maxnsplits, stop.rule, p.crit, det)

  if (det == "after" | det == "during") {
    y <- y - t(w[1:ncol(x)] %*% t(x))
    if (!is.null(exog)) w[-c(1:ncol(x))] <- coefficients(lm(y~ -1 + exog))
  }

  return(w)

}

