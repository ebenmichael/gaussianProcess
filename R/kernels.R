 ## Various kernel functions for use with GP regression


#' Compute the pairwise differences between X and Y
#' @param X An nx x d matrix
#' @param Y An ny x d matrix
#'
#' @return An nx * ny / 2 x d matrix of the pairwise distances in all dimensions
pdiff <- function(X, Y) {
    nx <- dim(X)[1]
    ny <- dim(Y)[1]
    # get the indices for the distance matrix
    pts <- expand.grid(1:nx, 1:ny)
    pts <- as.matrix(cbind(pts[,1], pts[,2]))
    # only keep upper triangle indices
    #pts <- matrix(pts[pts[,2] >= pts[,1],], ncol=2)
    return(X[pts[, 1], drop=F] - Y[pts[,2], drop=F])
}


#' Compute the scaled pairwise distances between X and Y
#' @param X An nx x d matrix
#' @param Y An ny x d matrix
#' @param scales The scale parameters
#'
#' @return An nx x ny matrix of scaled distances
pdist_scaled <- function(X, Y, scales) {
    # pairwise distances
    nx <- dim(X)[1]
    ny <- dim(Y)[1]
    dx <- dim(X)[2]
    dy <- dim(Y)[2]
    
    if(dx != dy) stop("Matrices have different dimensions") else d <- dx
    # get pairwise differences, square and scale
    diff.mat <- pdiff(X,Y) ^ 2
    diff.mat <- t(t(diff.mat) / scales^2)
    distances <- rowSums(diff.mat)
    #dist.mat <- matrix(0, nrow=nx, ncol=ny)
    #utri <- upper.tri(dist.mat, diag=TRUE)
    #dist.mat[utri] <- distances
    #ltri <- lower.tri(dist.mat, diag=TRUE)
    #dist.mat[ltri] <- distances
    dist.mat <- matrix(distances, nrow=dim(X)[1])
    return(dist.mat)
}


#' Compute the RBF (Squared Exponential) kernel between two matrices
#'
#' @param X An n x d matrix
#' @param Y An n x d matrix
#' @param scales The scale parameters, defaults to all ones
#' @param amplitude A parameter which is multiplied by the kernel, defaults to 1
#'
#' @return An n x n kernel matrix formed from \code{X} and \code{Y}
#'
#' @export
rbf <- function(X, Y, scales=NULL, amplitude=1){
    d <- dim(X)[1]
    if(is.null(scales)) {
        scales = rep(1, d)
    }
    distMat <- pdist_scaled(X, Y, scales)
    # Get the squared exponential kernel matrix
    K <- exp(- distMat / 2) * amplitude^2
    return(K)
    
}

#' Compute the Matern kernel between two matrices
#'
#' @param X An n x d matrix
#' @param Y An n x d matrix
#' @param scales The scale parameters, defaults to all ones
#' @param order whether to compute a 3/2 or 5/2 kernel, defaults to 5/2
#' @param amplitude A parameter which is multiplied by the kernel, defaults to 1
#' @return An n x n kernel matrix formed from \code{X} and \code{Y}
matern <- function(X, Y, scales=NULL, order=5/2, amplitude=1){
    # compute pairwise distances
    d <- dim(X)[1]
    if(is.null(scales)) {
        scales = rep(1, d)
    }
    distMat <- pdist_scaled(X, Y, scales)
    if(order == 3/2){
        K <- (1 + sqrt(3 * distMat)) * exp(- sqrt(3) * distMat)
        K <- K * amplitude^2
    }
    else if(order == 5/2) {
        K <- (1 + sqrt(5 * distMat) + 5/3 * distMat ) *
              exp(-sqrt(5 * distMat)) * amplitude^2
    }
    else{
        stop("Only Matern 3/2 and 5/2 are supported")
    }
    
    return(K)
}
