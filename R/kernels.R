 ## Various kernel functions for use with GP regression

#' Compute the RBF (Squared Exponential) kernel between two matrices
#'
#' @param X An n x d matrix
#' @param Y An n x d matrix
#' @param scale The scale parameter, defaults to 1
#' @param amplitude A parameter which is multiplied by the kernel, defaults to 1
#'
#' @return An n x n kernel matrix formed from \code{X} and \code{Y}
#'
#' @export
rbf <- function(X, Y, scale=1, amplitude=1){
    # pairwise distances
    n <- dim(X)[1]
    distMat <- matrix(proxy::dist(X, Y), nrow=dim(X)[1])
    # square to get the squared norm
    distMat <- distMat ^ 2
    # Get the squared exponential kernel matrix
    K <- exp(- distMat / (2 * scale ^ 2)) * amplitude
    return(K)
    
}

#' Compute the Matern kernel between two matrices
#'
#' @param X An n x d matrix
#' @param Y An n x d matrix
#' @param scale The scale parameter, defaults to 1
#' @param order whether to compute a 3/2 or 5/2 kernel, defaults to 5/2
#' @param amplitude A parameter which is multiplied by the kernel, defaults to 1
#' @return An n x n kernel matrix formed from \code{X} and \code{Y}
matern <- function(X, Y, scale=1, order=5/2, amplitude=1){
    # compute pairwise distances
    distMat <- matrix(proxy::dist(X, Y), nrow=dim(X)[1])
    if(order == 3/2){
        K <- (1 + sqrt(3) * distMat / scale) * exp(- sqrt(3) * distMat / scale)
        K <- K * amplitude
    }
    else if(order == 5/2) {
        K <- (1 + sqrt(5) * distMat / scale + (5 * distMat ^ 2) / (3 * scale ^ 2)) *
              exp(- (sqrt(5) * distMat) / scale) * amplitude
    }
    else{
        stop("Only Matern 3/2 and 5/2 are supported")
    }
    
    return(K)
}
