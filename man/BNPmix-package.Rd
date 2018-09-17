\name{BNPmix-package}
\alias{BNPmix-package}
\alias{BNPmix}
\docType{package}
\title{
\packageTitle{BNPmix}
}
\description{
\packageDescription{BNPmix}
Conditional method for Dirichlet process mixture models
}
\details{
\packageDESCRIPTION{BNPmix}
\packageIndices{BNPmix}
This package offers a fast C++ implementation to estimate via conditional methods Dirichlet process mixture models. The models considered are of the form:

\deqn{X_i | \theta_i \sim k(X_i; \theta_i) \\ \theta_i|p \sim p \\ p \sim DP(\alpha P_0)}

We considered the univariate case, where \eqn{\theta = (\mu, \sigma^2)} and

\deqn{k(X;\mu, \sigma^2) \sim N(\mu, \sigma^2)}

with a Normal-Inverse-Gamma base measure

\deqn{P_0 \sim N(\mu; m_0,k_0 \cdot \sigma^2 ) \times IG(\sigma^2; a_0, b_0)}

and the multivariate case, where \eqn{\theta = (\matbf{\mu}, \Sigma)} and

\deqn{P_0 \sim N(\mathbf{\mu}; m_0,k_0 \cdot \Sigma ) \times IW(\Sigma; S_0, \nu_0)}.

The function \strong{BNPmixUS} estimate an univariate DPM model with Gaussian kernel via conditional Polya urn scheme.
The function \strong{condSLI} estimate an univariate DPM model with Gaussian kernel via conditional Slice sampler.
The function \strong{BNPmixUSmv} estimate a multivariate DPM model with Gaussian kernel via conditional Polya urn scheme.
The function \strong{condSLImv} estimate a multivariate DPM model with Gaussian kernel via conditional Slice sampler.


}
\author{
\packageAuthor{BNPmix}

Maintainer: \packageMaintainer{BNPmix}
}
\keyword{Bayesian, Bayesian nonparametric, MCMC}
\examples{
data_toy <- c(rnorm(100, -3, 1), rnorm(100, 3, 1))
grid <- seq(-7, 7, length.out = 50)
est_model <- BNPmixUS(data = data_toy, grid = grid, niter = 5000,
                    nburn = 1000, napprox = 100)
}
