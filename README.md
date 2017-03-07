# gaussianProcess
R package for Gaussian Process regression with various kernels


## Usage
If `X` is a matrix of training covariates and `y` a vector of training targets then you create a `gaussianProcess` and automatically tune the hyper parameters with various options (see doc) with

```R
gp <- gaussianProcess(X,y,options)
```

To predict the output for a new data matrix `X.new` run

```R
pred <- predict(gp, X.new)
```

Here `pred` is a list with the posterior mean in `pred$mean` and the posterior covariance in `pred$covariance`. If in 1 dimension, you can plot the posterior process between `x.min` and `x.max` with

```R
plot(gp, x.min, x.max)
```

### Using a mean function
If you define a mean function `mean.f` that takes in a matrix and returns a vector you can use this as the mean function for the GP with

```R
gp <- gaussianProcess(X,y, meanFunc=mean.f)
```
For instance you can use `glmnet` to train a mean function to pass into the GP as follows. Fit with glmnet (e.g. with lambda=1), get the beta and create the mean function:

```R
fit <- glmnet(x, y, lambda=1)
beta <- fit$beta
mean.f <- function(x) x %*% beta
```

Or you can use the `pryr` package to do a partial application

```R
fit <- glmnet(x, y, s=1)
mean.f <- pryr::partial(predict, object=fit)
```

### Cross Validation
You can Leave One Out CV errors (true - predicted) for all data points with `loocv(gp)`
