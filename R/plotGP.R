## Function to plot the results from a GP regression


#' Plot the predictive distribution of new data. Only works with 1 dimension
#'
#' @param gp gaussianProcess object
#' @param x.min Lowest x value to plot
#' @param x.max Highest x value to plot
#' @param by Spacing for x values, defaults to 0.1
#' @param plot.points Whether to plot the data as well, defaults to False
#' @param plot.mean Whether to plot the mean of the process
#' @param plot.int Whether to plot an interval around the mean
#' @param plot.samples Whether to plot samples
#' @param plot.predictive Whether to plot the posterior predictive
#' @param int.size Posterior intervals are mean +- int.size * var, defaults to 1.96
#' @param n.samples Number of samples to plot
#'
#' @export
plot.gaussianProcess <- function(gp, x.min, x.max, by=0.1, 
                                 plot.points=TRUE, plot.mean=TRUE,
                                 plot.int=TRUE, plot.samples=FALSE,
                                 plot.predictive=FALSE, int.size=1.96,
                                 n.samples=10) {

    Xnew <- seq(x.min, x.max, by)
    # get the posterior predictive mean and covariance
    pred <- predict(gp, Xnew)

    # put the posterior predictive in a data frame
    vars <- diag(pred$covariance)
    if(plot.predictive) {
        vars <- vars + gp$noise.var
    }
    upper <- pred$mean + int.size * vars
    lower <- pred$mean - int.size * vars
    df.pred <- data.frame(x=Xnew, mean=pred$mean, upper=upper, lower=lower)
    if(plot.points) {
        # plot the original data
        df <- data.frame(x=gp$data, y=gp$target)
        gplt <- ggplot2::ggplot() +
            ggplot2::geom_point(data=df, ggplot2::aes(x=x, y=y), alpha=.5)
    }
    else {
        gplt <- ggplot2::ggplot(df.pred, ggplot2::aes(x=x, y=mean))
    }

    # plot the posterior mean
    if(plot.mean) {
        gplt <- gplt + ggplot2::geom_line(data=df.pred,
                                          ggplot2::aes(x=x,y=mean),
                                          color='red')
    }
    if(plot.int) {
        gplt <- gplt + ggplot2::geom_ribbon(data=df.pred,
                       ggplot2::aes(x=x, ymin=lower, ymax=upper),
                       alpha=.2)
    }

    # plot the samples if desired
    if(plot.samples) {
        samples <- data.frame(t(sample(gp, Xnew, n.samples)))
        samples$x <- Xnew
        samples <- reshape2::melt(samples, id="x")
        # make the plots more transparent if also plotting the mean
        if(plot.mean) alpha=.3 else alpha=1
        gplt <- gplt + ggplot2::geom_line(data=samples,
                                          ggplot2::aes(x=x, y=value,
                                                       color=variable),
                                          alpha=alpha) +
            ggplot2::theme(legend.position="none")
    }
    # show the plot
    gplt
    
}
