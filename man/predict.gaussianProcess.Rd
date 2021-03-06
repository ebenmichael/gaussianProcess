% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predictGP.R
\name{predict.gaussianProcess}
\alias{predict.gaussianProcess}
\title{Get the posterior distribution for f(X)}
\usage{
\method{predict}{gaussianProcess}(gp, X)
}
\arguments{
\item{gp}{A gaussianProcess object with data, kernel, etc.}

\item{X}{An n x d matrix of new observations}
}
\value{
posterior A list with fields mean and covariance
}
\description{
Get the posterior distribution for f(X)
}

