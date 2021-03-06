% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_generic.R
\name{BNPdens}
\alias{BNPdens}
\title{BNPdens class constructor}
\usage{
BNPdens(
  density = NULL,
  data = NULL,
  grideval = NULL,
  grid_x = NULL,
  grid_y = NULL,
  clust = NULL,
  mean = NULL,
  beta = NULL,
  sigma2 = NULL,
  probs = NULL,
  niter = NULL,
  nburn = NULL,
  tot_time = NULL,
  univariate = TRUE,
  regression = FALSE,
  dep = FALSE,
  group_log = NULL,
  group = NULL,
  wvals = NULL
)
}
\arguments{
\item{density}{a matrix containing the values taken by the density at the grid points;}

\item{data}{a dataset;}

\item{grideval}{a set of values where to evaluate the density;}

\item{grid_x}{regression grid, independent variable;}

\item{grid_y}{regression grid, dependent variable;}

\item{clust}{a (\code{niter - nburn}) \eqn{\times}{x} \code{nrow(data)}-dimensional matrix containing
the cluster labels for each observation (cols) and MCMC iteration (rows);}

\item{mean}{values for the location parameters;}

\item{beta}{coefficients for regression model (only for \code{PYregression});}

\item{sigma2}{values of the scale parameters;}

\item{probs}{values for the mixture weights;}

\item{niter}{number of MCMC iterations;}

\item{nburn}{number of MCMC iterations to discard as burn-in;}

\item{tot_time}{total execution time;}

\item{univariate}{logical, \code{TRUE} if the model is univariate;}

\item{regression}{logical, \code{TRUE}  for the output of  \code{PYregression};}

\item{dep}{logical, \code{TRUE}  for the output of  \code{DDPdensity};}

\item{group_log}{group allocation for each iteration (only for \code{DDPdensity});}

\item{group}{vector, allocation of observations to strata (only for \code{DDPdensity});}

\item{wvals}{values of the processes weights (only for \code{DDPdensity}).}
}
\description{
A constructor for the \code{BNPdens} class. The class \code{BNPdens} is a named list containing
the output generated by a specified Bayesian nonparametric mixture model implemented by means of
a specified MCMC strategy, as in \code{PYdensity},  \code{DDPdensity}, and  \code{PYregression}.
}
\examples{
data_toy <- c(rnorm(100, -3, 1), rnorm(100, 3, 1))
grid <- seq(-7, 7, length.out = 50)
est_model <- PYdensity(y = data_toy, mcmc = list(niter = 100,
                      nburn = 10, nupd = 100), output = list(grid = grid))
str(est_model)
class(est_model)
}
