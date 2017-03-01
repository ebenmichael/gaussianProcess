## Function to predict new data points for a GP


#' Get the posterior distribution for f(X)
#'
#' @param gp A gaussianProcess object with data, kernel, etc.
#' @param X An n x d matrix of new observations
#' @return posterior A list with fields mean and covariance
#'
#' @export
predict.gaussianProcess <- function(gp, X) {

    if(! is.matrix(X)){
        X <- as.matrix(X)
    }
    
    Knew <- gp$kernel(X, gp$data)
    # get the posterior mean
    means <- gp$meanFunc(X)
    postMean <- means + Knew %*% gp$alpha

    Kself <- gp$kernel(X, X)
    Linv <- solve(gp$chol)
    postCov <- Linv %*% t(Knew)
    postCov <- (Knew %*% t(Linv)) %*% postCov
    postCov <- Kself - postCov
    # add prior noise
    # postCov <- postCov + gp$noiseVar * diag(dim(postCov)[1])
    posterior <-  list(mean=postMean, covariance=postCov)
    return(posterior)
    
}

## Override the sample method to be able to have an S3 generic (funky stuff)
sample.default = sample
sample <- function(obj, ...) {
    UseMethod("sample")
}

#' Sample from the posterior distribution for f(X) (or the posterior predictive)
#'
#' @param gp A gaussianProcess object with data, kernel, etc.
#' @param X An n x d matrix of new observations
#' @param n Number of samples
#' @param pred Boolean, whether to sample from the posterior predictive
#' @return An n dimensional vector of samples from the posterior (predictive)
#'
#' @export
sample.gaussianProcess <- function(gp, X, n, pred=F) {
    # get the posterior of f
    post <- predict(gp, X)
    mu <- post$mean
    Sigma <- post$covariance
    if(pred) {
        Sigma <- Sigma + gp$noiseVar * diag(dim(Sigma)[1])
    }
    # sample
    y <- MASS::mvrnorm(n, mu, Sigma)
    return(y)
}
