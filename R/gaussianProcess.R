## Function to create a Gaussian Process object with a given mean function,
## kernel, and data


zeroFunction <- function(x) 0

#' Create a gaussianProcess object with a given mean function and covariance
#' kernel, and data
#' @param X An n x d matrix of covariates for observed data
#' @param y An n dimensional vector of outputs for observed data
#' @param meanFunc The mean function of the process, defaults to 0
#' @param kernel The covariance kernel of the process, defaults to rbf
#' @param hyper.params The hyper parameters for the kernel
#'                     (\code{c(amplitude, scales)}), if NULL then optimize them
#'                     using the log likelihood
#' @param noise.var The variance of the noise around the function
#' @param order The order of the kernel, defaults to 5/2
#' @param verbose Level of information printed out, defaults to 0
#'
#' @return A gaussianProcess with these options
#'
#' @export
gaussianProcess <- function(X, y, meanFunc=zeroFunction, kernel=rbf,
                            hyper.params=NULL, noise.var=1, order=5/2,
                            verbose=0) {
    if(! is.matrix(X)) {
        X <- as.matrix(X)
    }
    d <- dim(X)[2]
    # center the data by the mean function
    y.centered <- y - meanFunc(X)
    # optimize the hyper parameters if necessary
    if(is.null(hyper.params)) {
        hyper.params <- optimize_hyper_params(X, y.centered, kernel, order,
                                              noise.var, verbose=verbose)
    }
    n.params <- d + 1
    # initialize the list
    gaussianProcess <- list()
    # set the kernel
    gaussianProcess$kernel.type <- kernel
    gaussianProcess$kernel.order <- order
    gaussianProcess$hyper.params <- hyper.params
    if(identical(kernel, rbf)) {
        gaussianProcess$kernel <- pryr::partial(kernel,
                                                scales=hyper.params[2:n.params],
                                                amplitude=hyper.params[1])
    }
    else if(identical(kernel, matern)) {
        gaussianProcess$kernel <- pryr::partial(kernel,
                                                order=order,
                                                scales=hyper.params[2:n.params],
                                                amplitude=hyper.params[1])                                                
    }
    else {
        stop("Kernel must be rbf or matern. Working to fix this")
    }

    gaussianProcess$meanFunc <- meanFunc
    
    # get the kernel matrix and cholesky decomposition for future prediction
    K <- gaussianProcess$kernel(X, X)
    L <- t(chol(K + noise.var * diag(dim(X)[1])))
    alpha <- backsolve(t(L), forwardsolve(L, y.centered))
    # keep around all this stuff
    gaussianProcess$cholesky <- L
    gaussianProcess$alpha <- alpha
    gaussianProcess$data <- X
    gaussianProcess$target <- y
    gaussianProcess$noise.var <- noise.var
    class(gaussianProcess) <- "gaussianProcess"
    return(gaussianProcess)
}
