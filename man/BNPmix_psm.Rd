% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{BNPmix_psm}
\alias{BNPmix_psm}
\title{C++ function - compute the posterior similarity matrix}
\arguments{
\item{M}{a matrix (r x n), r number of replications, n number of observations}
}
\description{
C++ function - compute the posterior similarity matrix
}
\examples{
{
  M <- matrix(c(1,1,1,2,1,1,2,2,1,1,2,1,1,2,1,1), ncol = 4)
  BNPmix_psm(M)
}

}
\keyword{internal}
