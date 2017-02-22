# gaussianProcess
R package for Gaussian Process regression with various kernels


## Usage
If `X` is a matrix of training covariates and `y` a vector of training targets then you create a `gaussianProcess` with various options (see doc) with

```gp <- gaussianProcess(X,y,options)```

To predict the output for a new data matrix `X.new` run

```pred <- predict(gp, X.new)```

Here `pred` is a list with the posterior mean in `pred$mean` and the posterior covariance in `pred$covariance`. If in 1 dimension, you can plot the posterior process between `x.min` and `x.max` with

```plot(gp, x.min, x.max)```

