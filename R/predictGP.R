## Function to predict new data points for a GP

source("gaussianProcess.R")

#' Get the predictive distribution for new data X
#'
#' @param gp A gaussianProcess object with data, kernel, etc.
#' @param X An n x d matrix of new observations
#' @return posterior A list with fields mean and covariance
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
    posterior <-  list(mean=postMean, covariance=postCov)
    return(posterior)
    
}
