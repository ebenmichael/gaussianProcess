## Routines and helper functions to perform hyper parameter optimization


############Computing the objective and derivatives#############################


## Computing the log marginal likelihood
#' Compute the marginal log likelihood for a GP with given data, kernel, etc.
#'
#' @param hyper.params Hyper parameters as a vector (for rbf c(amplitude, scales))
#' @param X An n x d matrix of covariates for observed data
#' @param y An n dimensional vector of outputs for observed data
#' @param kernel.type The kernel function to use
#' @param order The order of the kernel, defaults to 5/2
#' @param noise.var The variance of the noise around the function
#'
#' @return Log marginal likelihood
gp_log_marg_like <- function(hyper.params, X, y, kernel.type=rbf,
                                  order=5/2, noise.var=1) {
    # apply hyper parameters to kernel
    n.params <- length(hyper.params)
    if(identical(kernel.type, rbf)) {
        kernel <- pryr::partial(kernel.type, scales=hyper.params[2:n.params],
                                   amplitude=hyper.params[1])
    }
    else if(identical(kernel.type, matern)) {
        kernel <- pryr::partial(kernel.type, scales=hyper.params[2:n.params],
                                order=order, amplitude=hyper.params[1])
    }
    # precompute the cholesky factors for the matrix inversion
    K <- kernel(X, X) + noise.var * diag(dim(X)[1])
    L <- t(chol(K))
    alpha <- backsolve(t(L), forwardsolve(L, y))
    # compute the log marginal likelihood
    return(gp_log_marg_like_helper(y, K, alpha))
}


#' Compute the log marginal likelihood from parameters of the GP
#'
#' @param y An n dimensional vector of outputs for observed data
#' @param K The n x n covariance matrix on the seen data
#' @param alpha An n dimensional vector, alpha = K^-1 y
#'
#' @return Log marginal likelihood
gp_log_marg_like_helper <- function(y, K, alpha) {
    # compute the parts of the log marginal likelihood
    n <- dim(K)[1]
    logp1 <- -1/2 * t(y) %*% alpha # data fit
    logp2 <- -1/2 * log(det(K)) # complexity penalty
    logp3 <- -n/2 * log(2 * pi) # normalization constant
    return(logp1 + logp2 + logp3)
}

## Generic S3 method
log_marginal_like <- function(obj, ...) {
    UseMethod("log_marginal_like")
}

#' Compute the log marginal likelihood for a GP
#' @param gp A gaussianProcess object
#'
#' @return Log marginal likelihood
log_marginal_like.gaussianProcess <- function(gp) {
    K <- gp$cholesky %*% t(gp$cholesky)
    y.centered <- gp$target - gp$meanFunc(gp$data)
    return(gp_log_marg_like_helper(y.centered, K, gp$alpha))
}


## Derivatives for optimization


#' Compute the derivative of the log marginal likelihood
#' @param hyper.params Hyper parameters as a vector (for rbf c(amplitude, scales))
#' @param X An n x d matrix of covariates for observed data
#' @param y An n dimensional vector of outputs for observed data
#' @param kernel.type The kernel function to use
#' @param order The order of the kernel, defaults to 5/2
#' @param noise.var The variance of the noise around the function
#'
#' @return The gradient of the log marginal likelihood at \code{hyper.params}
gp_lml_deriv <- function(hyper.params, X, y, kernel.type, order=5/2, noise.var=1) {
    n.params <- length(hyper.params)
    if(identical(kernel.type, rbf)) {
        kernel <- pryr::partial(kernel.type, scales=hyper.params[2:n.params],
                                amplitude=hyper.params[1])
        kernel.deriv <- rbf_deriv
    }
    else if(identical(kernel.type, matern)) {
        kernel <- pryr::partial(kernel.type, scales=hyper.params[1:n.params-1],
                                order=order, amplitude=hyper.params[n.params])
    }
    # precompute the cholesky factors for the matrix inversion
    K <- kernel(X, X) + noise.var * diag(dim(X)[1])
    L <- t(chol(K))
    alpha <- backsolve(t(L), forwardsolve(L, y))
    K.inv <- backsolve(t(L), forwardsolve(L, diag(dim(K)[1])))
    diff.mat <- pdiff(X, X)
    return(gp_lml_deriv_helper(hyper.params, diff.mat, y,
                               K, K.inv, alpha, kernel.deriv, noise.var))
}

                                                                     
#' Compute the derivative of the log marginal likelihood (helper function)
#' @param hyper.params Hyper parameters as a vector (for rbf c(amplitude, scales))
#' @param diff.mat The result of \code{pdiff(X,X)} for input data \code{X}
#' @param y The target vector
#' @param K The kernel matrix of \code{X}
#' @param K.inv The inverse of \code{K} (precomputed)
#' @param alpha An n dimensional vector, alpha = K^-1 y
#' @param kernel.deriv The derivative function for the kernel
#' @param noise.var The variance of the noise around the function
#'
#' @return The gradient of the log marginal likelihood at \code{hyper.params}
gp_lml_deriv_helper <- function(hyper.params, diff.mat, y,
                                K, K.inv, alpha, kernel.deriv, noise.var) {
    alpha.sq <- alpha %*% t(alpha)
    mat1 <- alpha.sq - K.inv
    d <- length(hyper.params)
    grad <- numeric(d)
    # compute gradient for each hyper parameter 
    # TODO: better than for loop?
    for(j in 1:d) {
        # compute jth partial derivative at kernel matrix
        dKdt <- kernel.deriv(j, hyper.params, diff.mat, K, noise.var)
        # only compute the diagonals of this matrix since we compute the trace
        grad[j] <- 0.5 * sum(diag((colSums(t(mat1) * dKdt))))
    }
    return(grad)
}



## Derivatives of kernel functions
#' Compute the jth partial derivative  of the rbf kernel 
#' @param j The coordinate of the derivative to compute
#' @param hyper.params Hyper parameters as a vector (for rbf c(amplitude, scales))
#' @param diff.mat The result of \code{pdiff(X,X)} for input data \code{X}
#' @param K The kernel matrix of \code{X}
#' @param noise.var The variance of the noise around the function
#'
#' @return The jth partial derivative at \code{hyper.params}
rbf_deriv <- function(j, hyper.params, diff.mat, K, noise.var) {
    n <- dim(K)[1]
    if(j == 1) {
        # derivative of amplitude
        return(2 * (K - noise.var * diag(n)) / hyper.params[1])
    } else {
        # compute the matrix of differences in jth coordinate
        # deal with the onde dimensional case
        if(is.null(dim(diff.mat))) {
            dist.mat <- matrix(diff.mat^2, nrow=n)
        } else {
            dist.mat <- matrix(diff.mat[,j-1]^2, nrow=n)
        }
        return((K - noise.var * diag(n)) *  dist.mat / hyper.params[j]^3)
    }
}

#############################Optimization#############################

#' Optimize the hyper parameters of a GP
#' @param X An n x d matrix of covariates for observed data
#' @param y An n dimensional vector of outputs for observed data
#' @param kernel The covariance kernel of the process, defaults to rbf
#' @param order The order of the kernel, defaults to 5/2
#' @param noise.var The variance of the noise around the function
#' @param verbose Level of information printed out, defaults to 0
#'
#' @return Hyper parameters which optimize the log marginal likelihood
optimize_hyper_params <- function(X, y, kernel.type=rbf,
                                  order=5/2, noise.var=1, verbose=0) {

    obj <- pryr::partial(gp_log_marg_like, X=X, y=y,
                         kernel.type=kernel.type,
                         order=order, noise.var=noise.var)
    grad <- pryr::partial(gp_lml_deriv, X=X, y=y, kernel.type=kernel.type,
                          order=order, noise.var=noise.var)

    n.params <- dim(X)[2] + 1
    # initialize at ones
    theta0 <-rep(1, n.params)
    # optimize
    ctrl <- list()
    ctrl$trace <- verbose
    ctrl$fnscale <- -1
    opt <- stats::optim(theta0, obj, grad, method="BFGS", control=ctrl)
    return(opt$par)
}

