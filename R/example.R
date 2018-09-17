# # test
#
# library(vars)
# library(HVAR)
#
# data(Canada)
# y <- apply(Canada, 2, function(x) (x-min(x))/(max(x)-min(x)))
#
# mod1 <- HVAR(y, type = "none", distance.lambda = 0, det = "during")
# mod2 <- HVAR2(y, type = "none", distance.lambda = 0)
# mod3 <- VAR(y, type = "none")
#
# Bcoef(mod1)
# Bcoef(mod2)
# Bcoef(mod3)
#
# y <- Canada
# p <- 1
# type <- "none"
# clust.method <- "AGNES"
# clust.control <- NULL
# thresh.val <- 1
# stop.rule <- 1
# maxnsplits <- 1
#
# xa <- rnorm(100)
# xb <- rnorm(100) + xa
# xc <- rnorm(100) - xb
# xd <- rnorm(100) + xa
# y <- data.frame(cbind(xa,xb,xc,xd))
# cor(y)
