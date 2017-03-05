## Function to create a Gaussian Process object with a given mean function,
## kernel, and data


zeroFunction <- function(x) 0

#' Create a gaussianProcess object with a given mean function and covariance
#' kernel, and data
#' @param X An n x d matrix of covariates for observed data
#' @param y An n dimensional vector of outputs for observed data
#' @param meanFunc The mean function of the process, defaults to 0
#' @param kernel The covariance kernel of the process, defaults to rbf
#' @param noiseVar The variance of the noise around the function
#' @param scales The scales of the kernel, defaults to ones
#' @param order The order of the kernel, defaults to 5/2
#' @param amplitude A parameter which is multiplied by the kernel, defaults to 1
#'
#' @return A gaussianProcess with these options
#'
#' @export
gaussianProcess <- function(X, y, meanFunc=zeroFunction, kernel=rbf,
                            noiseVar=1, scale=NULL, order=5/2, amplitude=1) {
    d <- dim(X)[1]
    if(is.null(scales)) {
        scales = rep(1, d)
    }
    # initialize the list
    gaussianProcess <- list()
    # set the kernel
    if(identical(kernel, rbf)) {
        gaussianProcess$kernel <- pryr::partial(kernel, scale=scale,
                                                amplitude=amplitude)
    }
    else if(identical(kernel, matern)) {
        gaussianProcess$kernel <- pryr::partial(kernel, scale=scale,
                                                order=order, amplitude=amplitude)
    }
    else {
        stop("Kernel must be rbf or matern. Working to fix this")
    }

    gaussianProcess$meanFunc <- meanFunc
    
    
    if(! is.matrix(X)) {
        X <- as.matrix(X)
    }
    # get the kernel matrix and cholesky decomposition for future prediction
    K <- gaussianProcess$kernel(X, X)
    L <- t(chol(K + noiseVar * diag(dim(X)[1])))
    # center the data by the mean function
    y.centered <- y - meanFunc(X)
    alpha <- solve(t(L), solve(L, y.centered))
    gaussianProcess$cholesky <- L
    gaussianProcess$alpha <- alpha
    gaussianProcess$data <- X
    gaussianProcess$target <- y
    gaussianProcess$noiseVar <- noiseVar
    class(gaussianProcess) <- "gaussianProcess"
    return(gaussianProcess)
}
