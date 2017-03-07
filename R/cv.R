## Routines for cross validation

## Generic S3 method
loocv <- function(object,...) {
    UseMethod("loocv")
}

#' Compute the leave one out CV estimate (Rasmussen & WIllaims pg 117)
#' @param gp The gaussianProcess
loocv.gaussianProcess <- function(gp) {
    # get the weight matrix on the observations
    L.inv <- forwardsolve(gp$cholesky, diag(dim(gp$cholesky)[1]))
    K.inv <- t(L.inv) %*% L.inv
    # compute the loocv values for each data point
    diffs <- gp$alpha / diag(K.inv)
    return(diffs)
}
