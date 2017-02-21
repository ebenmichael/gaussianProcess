## Function to plot the results from a GP regression

source("gaussianProcess.R")
source("predictGP.R")


#' Plot the predictive distribution of new data. Only works with 1 dimension
#'
#' @param gp gaussianProcess object
#' @param xnew Test data points to plot
#' @param int.size Posterior intervals are mean +- int.size * var, defaults to 1.96
#' @param y.var Whether to show predictive distribution for new y, defaults to False
#' @param plot.points Whether to plot the data as well, defaults to False
plot.gaussianProcess <- function(gp, Xnew, int.size=1.96, y.var=FALSE,
                                 plot.points=TRUE) {

    # get the posterior predictive mean and covariance
    pred <- predict(gp, Xnew)

    # put the posterior predictive in a data frame
    vars <- diag(pred$covariance)
    if(y.var) vars <- vars + gp$noiseVar
    upper <- pred$mean + int.size * vars
    lower <- pred$mean - int.size * vars
    df.pred <- data.frame(x=Xnew, mean=pred$mean, upper=upper, lower=lower)
    if(plot.points) {
        # plot the original data
        df <- data.frame(x=gp$data, y=gp$target)
        gplt <- ggplot2::ggplot() +
            ggplot2::geom_point(data=df, aes(x=x, y=y), alpha=.5)
    }
    else {
        gplt <- ggplot2::ggplot(df.pred, aes(x=x, y=mean))
    }

    # plot the posterior mean and the posterior intervals
    gplt <- gplt + geom_line(data=df.pred, aes(x=x,y=mean), color='red') +
        geom_ribbon(data=df.pred, aes(x=x, ymin=lower, ymax=upper), alpha=.3)
    # show the plot
    gplt
    
}