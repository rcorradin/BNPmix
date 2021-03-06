% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{MAR_mv_L}
\alias{MAR_mv_L}
\title{C++ function to estimate Pitman-Yor multivariate mixtures via marginal sampler - LOCATION}
\arguments{
\item{data}{a matrix of observations}

\item{grid}{matrix of points to evaluate the density}

\item{niter}{number of iterations}

\item{nburn}{number of burn-in iterations}

\item{m0}{expectation of location component}

\item{S20}{variance of location component}

\item{S0}{parameter of scale component}

\item{n0}{parameter of scale component}

\item{m1}{mean of hyperdistribtion of m0}

\item{k1}{scale factor of hyperdistribtion of m0}

\item{theta1}{df of hyperdistribtion of S20}

\item{Theta1}{matrix of hyperdistribution of S20}

\item{mass}{mass parameter}

\item{nupd}{number of iterations to show current updating}

\item{out_param}{if TRUE, return also the location and scale paramteres lists}

\item{out_dens}{if TRUE, return also the estimated density (default TRUE)}

\item{sigma_PY}{second parameter of PY}

\item{print_message}{print the status}

\item{light_dens}{if TRUE return only the posterior mean of the density}

\item{hyper, }{if TRUE use hyperpriors, default TRUE}
}
\description{
C++ function to estimate Pitman-Yor multivariate mixtures via marginal sampler - LOCATION
}
\keyword{internal}
