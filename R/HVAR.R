#' Estimation of a HVAR(p)
#'
#' Estimation of a HVAR using a clustered and diagonalised regressor matrix.
#'
#' This function implements a VAR model, where the parameter space is restricted along a hierarchical structure.
#'
#' The parameter 'clust.method' sets the cluster algorithm used to determine the hierarchical cluster. AGNES and DIANA implement the Agglomerative Nesting and Divisive Analysis algorithms of the 'cluster' package, while OPT using the package 'gaoptim' to optimise quasi-diagonality the regressor correlation matrix using the genetic permutation method 'GAPerm'.
#'
#' 'clust.control' allows additional arguments to be passed to the cluster method. Currently, only 2 are implemented: the AGNES method (default is single-linkage) and the maximum iterations (niter) for the genetic optimiser.
#'
#' The 'distance.lambda'controls the location where the cluster tree may be split into new branches. distance.lambda=1 searches for the optimal splitting location, while distance.lambda=0 permits any location. If max.splits>=ncol(x), where x is the regressor matrix, and distance.lambda=0, a normal VAR model is estimated.
#'
#' 'max.splits' controls the maximum splits to place at each level (default is 1). If there are more qualifying splits than max.splits, then the function prioritises tree symmetry.
#'
#' @param y Data item containing the endogenous variables
#' @param p Integer for the lag order (default is p=1)
#' @param type Type of deterministic regressors to include
#' @param det Should deterministic items be partialled out before, after or during coefficient estimation
#' @param season Inclusion of centered seasonal dummy variables (integer value of frequency)
#' @param exogen Inclusion of exogenous variables
#' @param clust.method Hierarchical clustering algorithm used to diagonalise the regressor correlation matrix
#' @param clust.control Optional arguments to pass to clust.method
#' @param distance.lambda Splitting threshold (must be within the range 0-1). See details for more information
#' @param max.splits Maximum number of splits to place at each hierarchical level
#' @param min.cluster Minimum number of items at each final node of the hierarchy
#' @param ridge.lambda Add regularization parameter lambda in the calculation of restricted coefficients.

#' @return A 'varest' object (see package 'vars'). This can be used with methods such as 'Bcoef', 'predict', 'fevd', etc. from the 'vars'-package. Note that this is not a perfect solution. For instance, bootstrapping irf intervals will not work and the 'lm' objects contained in 'object$varest' are not native. The purpose is to enable the user to process results, but native functions will be added in future.

#' @export

HVAR <- function(
  y,
  p = 1,
  type = c("const", "trend", "both", "none"),
  season = NULL,
  exogen = NULL,
  clust.method = c("AGNES", "DIANA", "OPT"),
  clust.control = NULL,
  distance.lambda = 0,
  max.splits = 1,
  min.cluster = 1,
  ridge.lambda = 0
)
{

  y <- as.matrix(y)
  if (any(is.na(y))) stop("\nNAs in y.\n")
  if (ncol(y) < 2) stop("The matrix 'y' should contain at least two variables. For univariate analysis consider ar() and arima() in package stats.\n")
  if (is.null(colnames(y))) {
    colnames(y) <- paste("y", 1:ncol(y), sep = "")
    warning(paste("No column names supplied in y, using:",
                  paste(colnames(y), collapse = ", "), ", instead.\n"))
  }
  colnames(y) <- make.names(colnames(y))

  type <- match.arg(type)
  clust.method <- match.arg(clust.method)

  obs <- dim(y)[1]
  K <- dim(y)[2]
  sample <- obs - p

  ylags <- embed(y, dimension = p + 1)[, -(1:K)]
  colnames(ylags) <- sub(" ","", apply(expand.grid(colnames(y), 1:p), 1, paste, collapse=".l"))
  yend <- y[-c(1:p), ]

  if (type == "const") {
    rhs <- matrix(rep(1, sample), ncol = 1)
    colnames(rhs) <- "const"
  }
  else if (type == "trend") {
    rhs <- matrix(seq(p + 1, length = sample), ncol = 1)
    colnames(rhs) <- "trend"
  }
  else if (type == "both") {
    rhs <- cbind(rep(1, sample), seq(p + 1, length = sample))
    colnames(rhs) <- c("const", "trend")
  }
  else if (type == "none") {
    rhs <- NULL
  }
  if (!(is.null(season))) {
    season <- abs(as.integer(season))
    dum <- (diag(season) - 1/season)[, -season]
    dums <- dum
    while (nrow(dums) < obs) {
      dums <- rbind(dums, dum)
    }
    dums <- dums[1:obs, ]
    colnames(dums) <- paste("sd", 1:ncol(dums), sep = "")
    rhs <- cbind(rhs, dums[-c(1:p), ])
  }
  if (!(is.null(exogen))) {
    exogen <- as.matrix(exogen)
    if (!identical(nrow(exogen), nrow(y))) {
      stop("\nDifferent row size of y and exogen.\n")
    }
    if (is.null(colnames(exogen))) {
      colnames(exogen) <- paste("exo", 1:ncol(exogen),
                                sep = "")
      warning(paste("No column names supplied in exogen, using:",
                    paste(colnames(exogen), collapse = ", "), ", instead.\n"))
    }
    colnames(exogen) <- make.names(colnames(exogen))
    tmp <- colnames(ylags)
    ylags <- cbind(ylags, exogen[-c(1:p), ])
    colnames(ylags) <- c(tmp, colnames(exogen))
  }
  datamat <- as.data.frame(ylags)
  colnames(datamat) <- colnames(ylags)

  cor.mat <- cor(datamat)
  cov.mat <- cov(datamat)
  ord <- get.clusters(cor.mat, clust.method = clust.method, clust.control = clust.control, absolute = T)
  dim.x <- ncol(datamat) + ifelse(!is.null(rhs), ncol(rhs), 0)
  coefs <- matrix(NA, dim.x, K)
  rownames(coefs) <- c(colnames(datamat), colnames(rhs))
  colnames(coefs) <- colnames(y)

  bipart <- recbipart(datamat, exog = rhs, cov.mat, cor.mat, ord$order, ord$height, distance.lambda,
            maxnsplits = max.splits, stop.rule = min.cluster)
  XX <- bipart$XX
  sign.all <- bipart$sign.all
  ymat <- yend[,1:K]
  x.det <- cbind(bipart$x.all, rhs)
  d <- diag(ncol(XX))
  diag(d) <- bipart$level.counter
  beta.raw <- solve(XX - ridge.lambda*d, tol = .Machine$double.eps^2) %*% t(x.det) %*% ymat
  beta.det <- beta.raw[-c(1:ncol(bipart$x.all)),]
  beta.raw <- beta.raw[1:ncol(bipart$x.all),]
  for (i in 1:K) {
    beta <- colSums(apply(sign.all, 1, function(x) x * t(beta.raw[,i])))
    if (length(dim(beta.det))==2) coefs[,i] <- c(beta, beta.det[,i]) else coefs[,i] <- c(beta, beta.det[i])
  }

  # Create a varest object
  equation <- list()
  for (i in 1:K) {
    form <- paste(colnames(y)[i],"~ -1 +",paste(colnames(datamat), collapse = "+"))
    if (!is.null(rhs)) form <- c(form, paste("+", paste(colnames(rhs), collapse = "+")))
    yend_i <- data.frame(yend[,i])
    colnames(yend_i) <- colnames(yend)[i]
    temp.data <- if(!is.null(rhs)) data.frame(cbind(yend_i, datamat, rhs)) else data.frame(cbind(yend_i, datamat))
    temp <- lm(formula = form, data = temp.data)
    temp.data <- if(!is.null(rhs)) data.frame(cbind(datamat, rhs)) else data.frame(cbind(datamat))
    temp$coefficients <- coefs[,i]
    temp$residuals <- as.numeric(yend[,i] - as.matrix(temp.data) %*% coefs[,i])
    temp$fitted.values = as.numeric(as.matrix(temp.data) %*% coefs[,i])
    equation[[colnames(y)[i]]] <- temp
  }
  results <- list(varresult = equation,
                  datamat = data.frame(cbind(yend, temp.data)),
                  y = y,
                  type = type,
                  p = p,
                  K = K,
                  obs = sample,
                  totobs = obs,
                  restrictions = NULL,
                  call = match.call(),
                  bipart = bipart
  )
  class(results) <- "varest"
  return(results)

}
