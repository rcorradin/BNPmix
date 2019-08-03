#' @name PYdensity
#' @export PYdensity
#'
#' @title MCMC for Pitman-Yor mixtures of Gaussians
#' @description The \code{PYdensity} function generates a posterior density sample for a selection of univariate and multivariate Pitman-Yor
#' process mixture models with Gaussian kernels. See details below for the description of the different specifications of the implemented models.
#'
#'
#' @param y a  vector or matrix giving the data based on which the density is to be estimated;
#'
#' @param output a list of arguments for generating posterior output. It contains:
#' \itemize{
#' \item \code{grid}, a grid of points at which to evaluate the estimated posterior mean density; a data frame
#' obtained with the \code{expand.grid} function.
#'
#' \item \code{out_param}, if equal to \code{TRUE}, the function returns the draws of the kernel's
#'   parameters for each MCMC iteration, default is \code{FALSE}. See \code{value} for details.
#'
#' \item \code{out_type}, if \code{out_type = "FULL"}, the function returns the visited
#'  partitions and the realizations of the posterior density for each iterations.
#' If \code{out_type = "MEAN"}, the function returns the estimated partitions and the mean of the densities sampled at each iterations.
#' If \code{out_type = "CLUST"}, the function returns the estimated partition.
#' Default \code{"FULL"}.
#'
# \item \code{mcmc_dens}, if equal to \code{TRUE}, the function returns a total of \code{niter}-\code{nburn} realizations of the posterior
# density, that is one per stored iteration, evaluated at
# \code{grid} (default is \code{TRUE}). See \code{value} for details.
#
# \item \code{mean_dens}, if equal to \code{TRUE}, the function returns only the posterior mean of the density, default is \code{TRUE}. See \code{value} for details.
#' }
#'
#' @param mcmc a list of MCMC arguments:
#' \itemize{
#'   \item \code{niter} (mandatory), number of iterations.
#'
#'   \item \code{nburn} (mandatory), number of iterations to discard as burn-in.
#'
#'   \item \code{method}, the MCMC sampling method to be used (default is \code{'ICS'}). See \code{details}.
#'
#'   \item \code{model}, the type of model to be fitted (default is \code{'LS'}). See \code{details}.
#'
#'   \item \code{nupd}, argument controlling the number of iterations to be displayed on screen: the function reports
#'   on standard output every time \code{nupd} new iterations have been carried out (default is \code{niter/10}).
#'
#'   \item \code{print_message}, control option. If equal to \code{TRUE}, the status is printed
#'   to standard output every \code{nupd} iterations (default is \code{TRUE}).
#'
#'   \item \code{m_imp}, number of generated values for the importance sampling step of \code{method = 'ICS'} (default is 10). See \code{details}.
#'
#'   \item \code{hyper}, if equal to \code{TRUE}, hyperprior distributions on the base measure's
#'   parameters are added, as specified in \code{prior} and explained in \code{details} (default is \code{TRUE}).
#' }
#'
#'
#' @param prior a list giving the prior information. The list includes \code{strength} and \code{discount},
#' the strength and discount parameters of the Pitman-Yor process
#' (default are \code{strength = 1} and \code{discount = 0}, the latter leading to the Dirichlet process).
#' The remaining parameters depend on the model choice.
#' \itemize{
#'
#' \item If \code{model = 'L'} (location mixture) and \code{y} is univariate:
#'
#' \code{m0} and \code{s20} are
#'  mean and variance of the base measure on the location parameter (default are 0 and 1);
#'  \code{a0} and \code{b0} are shape and scale parameters of the inverse gamma prior on the common scale parameter
#'  (default are 2 and 1).
#'  If \code{hyper = TRUE},  optional hyperpriors on the base measure's parameters are added:
#'  specifically, \code{m1} and \code{k1} are the mean parameter and the scale factor defining the
#'  normal hyperprior on \code{m0} (default are 0 and 1), and
#'  \code{a1} and \code{b1} are shape and rate parameters of the inverse gamma hyperprior on \code{b0}
#'  (default are 1 and 2). See \code{details}.
#'
#' \item If \code{model = 'LS'} (location-scale mixture) and \code{y} is univariate:
#'
#'   \code{m0} and \code{k0} are the mean parameter and the scale factor defining the normal base measure
#'   on the location parameter (default are 0 and 1), and \code{a0} and \code{b0} are
#'   shape and scale parameters of the inverse gamma base measure on the scale parameter (default are 2 and 1).
#'  If \code{hyper = TRUE},  optional hyperpriors on the base measure's parameters are added:
#'  specifically, \code{m1} and \code{s21} are mean and variance parameters of the normal hyperprior on
#'  \code{m0} (default are 0 and 1);
#'  \code{tau1} and \code{zeta1} are shape and rate parameters of the gamma hyperprior on
#'  \code{k0} (default is 1 for both);
#'  \code{a1} and \code{b1} are shape and rate parameters of the gamma hyperprior on
#'  \code{b0}  (default is 1 for both). See \code{details}.
#'
#' \item If \code{model = 'L'} (location mixture) and \code{y} is multivariate (\code{p}-variate):
#'
#' \code{m0} and \code{S20} are
#'  mean and covariance of the base measure on the location parameter (default are a vector of zeroes and the identity matrix);
#'  \code{Sigma0} and \code{n0} are the parameters of the inverse Whishart prior on
#'   the common scale matrix.
#'  If \code{hyper = TRUE}, optional hyperpriors on the base measure's parameters are added:
#'  specifically, \code{m1} and \code{k1} are the mean parameter and the scale factor defining the
#'  normal hyperprior on \code{m0} (default are a vector of zeroes and 1), and
#'  \code{lambda} and \code{Lambda1} are the parameters (degrees of freedom and scale) of the inverse Wishart prior on \code{S20} (default are \code{p}+2 and the identity matrix). See \code{details}.
#'
#' \item If \code{model = 'LS'} (location-scale mixture) and \code{y} is multivariate (\code{p}-variate):
#'
#' \code{m0} and \code{k0} are the mean parameter and the scale factor defining the normal base measure on the
#' location parameter (default are a vector of zeros and 1), and
#' \code{n0} and \code{Sigma0} are the parameters (degrees of freedom and scale) of the inverse Wishart base measure on the location parameter
#' (default are \code{p}+2 and the identity matrix).
#'  If \code{hyper = TRUE}, optional hyperpriors on the base measure's parameters are added:
#'  specifically, \code{m1} and \code{S1} are mean and covariance matrix parameters of the normal hyperprior on
#'  \code{m0} (default are a vector of zeroes and the identity matrix);
#'  \code{tau1} and \code{zeta1} are shape and rate parameters of the gamma hyperprior on
#'  \code{k0} (default is 1 for both);
#'  \code{n1} and \code{Sigma1} are the parameters (degrees of freedom and scale) of the Wishart prior for \code{Sigma0} (default are \code{p}+2 and the identity matrix divided \code{p}+2). See \code{details}.
#'
#'
#' \item If \code{model = 'DLS'} (diagonal location-scale mixture):
#'
#' \code{m0} and \code{k0} are the mean vector parameter and the vector of scale factors defining the normal base measure
#' on the location parameter (default are a vector of zeroes and a vector of ones),
#' and \code{a0} and \code{b0} are vectors of
#' shape and scale parameters defining the base measure on the scale parameters (default are a vector of twos and a vector of ones).
#' If \code{hyper = TRUE}, optional hyperpriors on the base measure's parameters are added:
#'  specifically, \code{m1} and \code{s21} are vectors of mean and variance parameters for the normal hyperpriors on the components of
#'  \code{m0} (default are a vector of zeros and a vector of ones);
#'  \code{tau1} and \code{zeta1} are vectors of shape and rate parameters of the gamma hyperpriors on the components of
#'  \code{k0} (default is a vector of ones for both);
#'  \code{a1} and \code{b1} are vectors of shape and rate parameters of the gamma hyperpriors on the components of
#'  \code{b0}  (default is a vector of ones for both). See \code{details}.
#' }
#'
#'@details
#' This generic function fits a Pitman-Yor process mixture model for density estimation and clustering. The general model is
#' \deqn{\tilde f(y) = \int K(y; \theta) \tilde p (d \theta),} where \eqn{K(y; \theta)} is a kernel density with parameter
#' \eqn{\theta\in\Theta}. Univariate and multivariate Gaussian kernels are implemented with different specifications for the parametric space
#' \eqn{\Theta}, as described below.
#' The mixing measure \eqn{\tilde p} has a Pitman-Yor process prior with strength parameter \eqn{\vartheta},
#' discount parameter \eqn{\alpha}, and base measure \eqn{P_0} admitting the specifications presented below. For posterior sampling,
#' three MCMC approaches are implemented. See details below.
#'
#' \strong{Univariate data}
#'
#' For univariate \eqn{y} the function implements both a location and location-scale mixture model. The former assumes
#' \deqn{\tilde f(y) = \int \phi(y; \mu, \sigma^2) \tilde p (d \mu) \pi(\sigma^2),} where
#' \eqn{\phi(y; \mu, \sigma^2)} is a univariate Gaussian kernel function with mean \eqn{\mu} and variance \eqn{\sigma^2},
#' and \eqn{\pi(\sigma^2)} is an inverse gamma prior. The base measure is specified as
#' \deqn{P_0(d \mu) = N(d \mu; m_0, \sigma^2_0),}{%
#' P_0(d \mu) = N(d \mu; m0, s20),}
#'  and \eqn{\sigma^2 \sim IGa(a_0, b_0)}{\sigma^2 ~ IGa(a0, b0)}.
#'  Optional hyperpriors for the base measure's parameters are
#' \deqn{(m_0,\sigma^2_0) \sim N(m_1, \sigma^2_0 / k_1) \times IGa(a_1, b_1).}{%
#'       (m0,s20) ~ N(m1, s20 / k_1)  IGa(a1, b1).}
#'
#' The location-scale mixture model, instead, assumes
#' \deqn{\tilde f(y) = \int \phi(y; \mu, \sigma^2) \tilde p (d \mu, d \sigma^2)} with normal-inverse gamma base measure
#' \deqn{P_0 (d \mu, d \sigma^2) = N(d \mu; m_0, \sigma^2 / k_0) \times IGa(d \sigma^2; a_0, b_0).}{%
#'       P_0 (d \mu, d \sigma^2) = N(d \mu; m0, \sigma^2 / k0)  IGa(d \sigma^2; a0, b0)} and (optional) hyperpriors
#' \deqn{m_0 \sim N(m_1, \sigma_1^2 ),\quad k_0 \sim Ga(\tau_1, \zeta_1),\quad b_0 \sim Ga(a_1, b_1).}{%
#' m0 ~ N(m1, \sigma12 ),  k0 ~ Ga(\tau1, \zeta2),  b0 ~ Ga(a1, b1).}
#'
#'
#' \strong{Multivariate data}
#'
#' For multivariate \eqn{y} (\eqn{p}-variate) the function implements a location mixture model (with full covariance matrix) and two
#' different location-scale mixture models, with either full or diagonal covariance matrix. The location mixture model assumes
#' \deqn{\tilde f(y) = \int \phi_p(y; \mu, \Sigma) \tilde p (d \mu) \pi(\Sigma)} where
#' \eqn{\phi_p(y; \mu, \Sigma)} is a \eqn{p}-dimensional Gaussian kernel function with mean vector \eqn{\mu} and covariance matrix
#' \eqn{\Sigma}. The prior on \eqn{\Sigma} is inverse Whishart with parameters \eqn{\Sigma_0} and \eqn{\nu_0}, while the
#' base measure is
#' \deqn{P_0(d \mu) = N(d \mu;  m_0, S_0),}{P_0 (d \mu) = N(d \mu;  m0, S0),}
#'  with optional hyperpriors
#' \deqn{m_0 \sim N(m_1, S_0 / k_1),\quad S_0 \sim IW(\lambda_1, \Lambda_1).}{%
#'       m0 ~ N(m1, S0 / k1),  S0 ~ IW(\lambda1, \Lambda1).}
#'
#' The location-scale mixture model assumes
#'
#' \deqn{\tilde f(x) = \int \phi_p(y;  \mu, \Sigma) \tilde p (d  \mu, d \Sigma).} Two possible structures for \eqn{\Sigma}
#' are implemented, namely full and diagonal covariance. For the full covariance mixture model, the base measure is
#' the normal-inverse Wishart
#'  \deqn{P_0 (d \mu, d \Sigma) = N(d \mu;  m_0, \Sigma / k_0) \times IW(d \Sigma; \nu_0, \Sigma_0),}{%
#'        P_0 (d \mu, d \Sigma) = N(d \mu;  m0, \Sigma / k0)  IW(d \Sigma; \nu, \Sigma0),}
#'  with optional hyperpriors
#' \deqn{m_0 \sim N(m_1, S_1),\quad k_0 \sim Ga(\tau_1, \zeta_1),\quad b_0 \sim W(\nu_1, \Sigma_1).}{%
#'       m_0 ~ N(m1, S12), k0 ~ Ga(\tau1, \zeta1), b_0 ~ W(\nu1, \Sigma1).}
#'The second location-scale mixture model assumes a diagonal covariance structure. This is equivalent to write the
#'mixture model as a mixture of products of univariate normal kernels, i.e.
#' \deqn{\tilde f(y) = \int \prod_{r=1}^p \phi(y_r;  \mu_r, \sigma^2_r) \tilde p (d  \mu_1,\ldots,d \mu_p, d \sigma_1^2,\ldots,d \sigma_p^2).}
#' For this specification, the base measure is assumed defined as the product of \eqn{p} independent normal-inverse gamma distributions, that is
#'  \deqn{P_0 = \prod_{r=1}^p P_{0r}} where
#'  \deqn{P_{0r}(d \mu_r,d \sigma_r^2) = N(d \mu_r; m_{0r}, \sigma^2_r / k_{0r}) \times Ga(d \sigma^2_r; a_{0r}, b_{0r}).}{%
#'  P_{0r}(d \mu_r, d \sigma_r^2) = N(d \mu_r; m_{0j}, \sigma^2_r / k_{0r}) Ga(d \sigma^2_r; a_{0r}, b_{0r}). }
#'  Optional hyperpriors can be added, and, for each component, correspond to the set of hyperpriors considered
#'  for the univariate location-scale mixture model.
#'
#' \strong{Posterior simulation methods}
#'
#' This generic function implements three types of MCMC algorithms for posterior simulation.
#' The default method is the importance conditional sampler (Canale et al. 2019). Other options are
#' the marginal sampler (Neal, 2000)  and the dependent slice-efficient sampler (Kalli et al. 2011).
#' The importance conditional sampler performs an importance sampling step when updating the values of
#' individual parameters \eqn{\theta}, which requires to sample \code{m_imp} values from a suitable
#' proposal. Large values of \code{m_imp} are known to improve the mixing of the chain
#' at the cost of increased running time (Canale et al. 2019).
#'
#'
#' @return A \code{BNPdens} class object containing the estimated density and
#' the cluster allocations for each iterations. If \code{out_param = TRUE} the output
#' contains also the kernel specific parameters for each iteration. If \code{mcmc_dens = TRUE} the output
#' contains also a realization from the posterior density for each iteration. IF \code{mean_dens = TRUE}
#' the output contains just the mean of the realizations from the posterior density. The output contains
#' also informations as the number of iterations, the number of burn-in iterations, the used
#' computational time and the type of estimated model (\code{univariate = TRUE} or \code{FALSE}).
#'
#'
#' @examples
#' data_toy <- cbind(c(rnorm(100, -3, 1), rnorm(100, 3, 1)),
#'                   c(rnorm(100, -3, 1), rnorm(100, 3, 1)))
#' grid <- expand.grid(seq(-7, 7, length.out = 50),
#'                     seq(-7, 7, length.out = 50))
#' est_model <- PYdensity(y = data_toy, mcmc = list(niter = 200, nburn = 100),
#' output = list(grid = grid))
#' summary(est_model)
#' plot(est_model)
#'
#' @references
#'
#' Canale, A., Corradin, R., Nipoti, B. (2019), Importance conditional sampling for Bayesian nonparametric mixtures,
#' arXiv preprint, 	arXiv:1906.08147
#'
#' Kalli, M., Griffin,  J. E., and Walker,  S. G. (2011), Slice sampling mixture models.
#' Statistics and Computing 21, 93-105.
#'
#' Neal, R. M. (2000), Markov Chain Sampling Methods for Dirichlet Process Mixture Models,
#' Journal of Computational and Graphical Statistics 9, 249-265.
#'

PYdensity <- function(y,
                       mcmc = list(),
                       prior = list(),
                       output = list()){

  if(is.vector(y)){

    # mandatory parameters
    if(is.null(mcmc$niter)) stop("Missing number of iterations")
    if(is.null(mcmc$nburn)) stop("Missing number of burn-in iterations")

    # if mcmc misses some parts, add default
    niter = mcmc$niter
    nburn = mcmc$nburn
    method = ifelse(is.null(mcmc$method), "ICS", mcmc$method)
    model = ifelse(is.null(mcmc$model), "LS", mcmc$model)
    nupd = ifelse(is.null(mcmc$nupd), round(niter / 10), mcmc$nupd)
    m_imp = ifelse(is.null(mcmc$m_imp), 10, mcmc$m_imp)
    print_message = ifelse(is.null(mcmc$print_message), TRUE, mcmc$print_message)

    # output
    out_param = ifelse(is.null(output$out_param), FALSE, output$out_param)
    output$out_type = ifelse(is.null(output$out_type), "FULL", output$out_type)
    if(output$out_type == "FULL"){
      mcmc_dens = TRUE
      if(is.null(output$grid)){
        grid_use = 0
      } else {
        if(length(dim(output$grid)) > 1) stop("Wrong grid dimension")
        grid_use <- as.vector(output$grid)
      }
    } else if (output$out_type == "CLUST"){
      mcmc_dens = FALSE
      grid_use = 0
    }
    #
    # if(is.null(output$grid)){
    #   mcmc_dens = FALSE
    #   grid = 0
    # } else {
    #   # check the grid
    #   if(length(dim(output$grid)) > 1) stop("Wrong grid dimension")
    #   grid_use <- as.vector(output$grid)
    # }

    if(!(model == "LS" | model == "L")) stop("Wrong model setting")
    if(!(method == "ICS" | method == "SLI" | method == "MAR")) stop("Wrong method setting")

    hyper = ifelse(is.null(mcmc$hyper), TRUE, mcmc$hyper)
    # Check for different model parameters
    # if null, initialize default parameters
    if(model == "LS"){
      if(hyper){
        if(is.null(prior$a0)){ a0 = 2 } else { a0 = prior$a0 }
        if(is.null(prior$m1)){ m1 = 0 } else { m1 = prior$m1 }
        if(is.null(prior$s21)){ s21 = 1 } else { s21 = prior$s21 }
        if(is.null(prior$tau1)){ tau1 = 1 } else { tau1 = prior$tau1 }
        if(is.null(prior$zeta1)){ zeta1 = 1 } else { zeta1 = prior$zeta1 }
        if(is.null(prior$a1)){ a1 = 1 } else { a1 = prior$a1 }
        if(is.null(prior$b1)){ b1 = 1 } else { b1 = prior$b1 }

        if(is.null(prior$m0)){ m0 = rnorm(1, m1, sqrt(s21)) } else { m0 = prior$m0 }
        if(is.null(prior$k0)){ k0 = rgamma(1, tau1, zeta1) } else { k0 = prior$k0 }
        if(is.null(prior$b0)){ b0 = rgamma(1, a1, b1) } else { b0 = prior$b0 }
      } else {
        if(is.null(prior$m0)){ m0 = 0 } else { m0 = prior$m0 }
        if(is.null(prior$k0)){ k0 = 1 } else { k0 = prior$k0 }
        if(is.null(prior$a0)){ a0 = 2 } else { a0 = prior$a0 }
        if(is.null(prior$b0)){ b0 = 1 } else { b0 = prior$b0 }

        m1 <- s21 <- tau1 <- zeta1 <- a1 <- b1 <- 0
      }
    } else if(model == "L"){
      if(hyper){
        if(is.null(prior$a0)){ a0 = 2 } else { a0 = prior$a0 }
        if(is.null(prior$b0)){ b0 = 1 } else { b0 = prior$b0 }
        if(is.null(prior$m1)){ m1 = 0 } else { m1 = prior$m1 }
        if(is.null(prior$k1)){ k1 = 1 } else { k1 = prior$k1 }
        if(is.null(prior$a1)){ a1 = 2 } else { a1 = prior$a1 }
        if(is.null(prior$b1)){ b1 = 1 } else { b1 = prior$b1 }

        if(is.null(prior$m0)){ m0 = rnorm(1, m1, sqrt(b1 / ((a1 - 1) * k1))) } else { m0 = prior$m0 }
        if(is.null(prior$s20)){ s20 = 1 / rgamma(1, a1, 1/b1) } else { s20 = prior$s20 }
      } else {
        if(is.null(prior$m0)){ m0 = 0 } else { m0 = prior$m0 }
        if(is.null(prior$s20)){ s20 = 1 } else { s20 = prior$s20 }
        if(is.null(prior$a0)){ a0 = 2 } else { a0 = prior$a0 }
        if(is.null(prior$b0)){ b0 = 1 } else { b0 = prior$b0 }

        m1 <- k1 <- a1 <- b1 <- 0
      }
    }

    # process parameters
    strength = ifelse(is.null(prior$strength), 1, prior$strength)
    discount = ifelse(is.null(prior$discount), 0, prior$discount)
    if(strength < - discount) stop("strength must be greater than -discount")

    # estimate the model
    if(method == "ICS"){

      # call the ICS univariate function
      if(model == "LS"){
        est_model <- cICS(y, grid_use, niter, nburn, m0, k0, a0, b0, m1, s21, tau1, zeta1, a1, b1, strength,
                          m_imp, nupd, out_param, mcmc_dens, discount, print_message, hyper)
      } else if(model == "L"){
        est_model <- cICS_L(y, grid_use, niter, nburn, m0, s20, a0, b0, m1, k1, a1, b1, strength,
                            m_imp, nupd, out_param, mcmc_dens, discount, print_message, hyper)
      }
    } else if(method == "SLI"){

      # call the SLI univariate function
      if(model == "LS"){
        est_model <- cSLI(y, grid_use, niter, nburn, m0, k0, a0, b0, m1, s21, tau1, zeta1, a1, b1,
                          strength, nupd, out_param, mcmc_dens, discount, print_message, hyper)
      } else if(model == "L"){
        est_model <- cSLI_L(y, grid_use, niter, nburn, m0, s20, a0, b0, m1, k1, a1, b1,
                            strength, nupd, out_param, mcmc_dens, discount, print_message, hyper)
      }
    } else if(method == "MAR"){

      # call the MAR univariate function
      if(model == "LS"){
        est_model <- MAR(y, grid_use, niter, nburn, m0, k0, a0, b0, m1, s21, tau1, zeta1, a1, b1,
                         strength, nupd, out_param, mcmc_dens, discount, print_message, hyper)
      } else if(model == "L"){
        est_model <- MAR_L(y, grid_use, niter, nburn, m0, s20, a0, b0, m1, k1, a1, b1,
                           strength, nupd, out_param, mcmc_dens, discount, print_message, hyper)
      }

    }

    # return the results
    if(!isTRUE(out_param)){
      if(isTRUE(mcmc_dens)){
        result <- BNPdens(density = est_model$dens,
                          data = y,
                          grideval = grid_use,
                          clust = est_model$clust,
                          niter = niter,
                          nburn = nburn,
                          tot_time = est_model$time,
                          univariate = TRUE)
      }else{
        result <- BNPdens(data = y,
                          clust = est_model$clust,
                          niter = niter,
                          nburn = nburn,
                          tot_time = est_model$time,
                          univariate = TRUE)
      }
    } else {
      if(isTRUE(mcmc_dens)){
        result <- BNPdens(density = est_model$dens,
                          data = y,
                          grideval = grid_use,
                          clust = est_model$clust,
                          mean = est_model$mu,
                          sigma2 = est_model$s2,
                          probs = est_model$probs,
                          niter = niter,
                          nburn = nburn,
                          tot_time = est_model$time,
                          univariate = TRUE)
      }else{
        result <- BNPdens(data = y,
                          clust = est_model$clust,
                          mean = est_model$mu,
                          sigma2 = est_model$s2,
                          probs = est_model$probs,
                          niter = niter,
                          nburn = nburn,
                          tot_time = est_model$time,
                          univariate = TRUE)
      }
    }
  } else if(!is.vector(y)){

    # Check the mandatory parameters
    if(is.null(mcmc$niter)) stop("Missing number of iterations")
    if(is.null(mcmc$nburn)) stop("Missing number of burn-in iterations")

    # if mcmc misses some parts, add default
    niter = mcmc$niter
    nburn = mcmc$nburn
    method = ifelse(is.null(mcmc$method), "ICS", mcmc$method)
    model = ifelse(is.null(mcmc$model), "LS", mcmc$model)
    nupd = ifelse(is.null(mcmc$nupd), round(niter / 10), mcmc$nupd)
    print_message = ifelse(is.null(mcmc$print_message), TRUE, mcmc$print_message)
    m_imp = ifelse(is.null(mcmc$m_imp), 10, mcmc$m_imp)

    # output
    out_param = ifelse(is.null(output$out_param), FALSE, output$out_param)
    output$out_type = ifelse(is.null(output$out_type), "FULL", output$out_type)
    if(output$out_type == "FULL"){
      mean_dens = FALSE
      mcmc_dens = TRUE
      if(is.null(output$grid)){
        mcmc_dens = FALSE
        grid_use = matrix(0, ncol = ncol(y))
      } else {
        # check the grid
        if(ncol(output$grid) != ncol(y)) stop("The dimensions of grid and data not match")
        grid_use <- as.matrix(output$grid)
      }
    } else if (output$out_type == "MEAN"){
      mean_dens = TRUE
      mcmc_dens = TRUE
      if(is.null(output$grid)){
        mcmc_dens = FALSE
        grid_use = matrix(0, ncol = ncol(y))
      } else {
        # check the grid
        if(ncol(output$grid) != ncol(y)) stop("The dimensions of grid and data not match")
        grid_use <- as.matrix(output$grid)
      }
    } else if (output$out_type == "CLUST"){
      mean_dens = FALSE
      mcmc_dens = FALSE
      grid_use = matrix(0, ncol = ncol(y))
    }

    if(!(model == "LS" | model == "L" | model == "DLS")) stop("Wrong model setting")
    if(!(method == "ICS" | method == "SLI" | method == "MAR")) stop("Wrong method setting")

    hyper = ifelse(is.null(mcmc$hyper), TRUE, mcmc$hyper)
    # Check for different model parameters
    if(model == "LS"){
      if(hyper){
        if(is.null(prior$m1)){ m1 = rep(0, ncol(y)) } else { m1 = prior$m1 }
        if(is.null(prior$S1)){ S1 = diag(1, ncol(y)) } else { S1 = prior$S1 }
        if(is.null(prior$tau1)){ tau1 = 1 } else { tau1 = prior$tau1 }
        if(is.null(prior$zeta1)){ zeta1 = 1 } else { zeta1 = prior$zeta1 }
        if(is.null(prior$n1)){ n1 = ncol(y) + 2 } else { n1 = prior$n1 }
        if(is.null(prior$Sigma1)){ Sigma1 = diag(1/n1, ncol(y)) } else { Sigma1 = prior$Sigma1 }
        if(is.null(prior$n0)){ n0 = ncol(y) + 2 } else { n0 = prior$n0 }

        if(is.null(prior$Sigma0)){ Sigma0 = rWishart(n = 1, Sigma = Sigma1, df = n1)[,,1] } else {Sigma0 = prior$Sigma0}
        if(is.null(prior$m0)){ m0 = as.vector(rnorm(ncol(y)) %*% t(chol(S1)) + m1)} else { m0 = prior$m0 }
        if(is.null(prior$k0)){ k0 = rgamma(1, tau1, 1 / zeta1)} else { k0 = prior$k0 }
      } else {
        if(is.null(prior$m0)){ m0 = rep(0, ncol(y)) } else { m0 = prior$m0 }
        if(is.null(prior$k0)){ k0 = 1 } else { k0 = prior$k0 }
        if(is.null(prior$Sigma0)){ Sigma0 = diag(1, ncol(y)) } else { Sigma0 = prior$Sigma0 }
        if(is.null(prior$n0)){ n0 = ncol(y) + 2 } else { n0 = prior$n0 }

        m1 <- rep(0, ncol(y))
        n1 <- tau1 <- zeta1 <- 0
        Sigma1 <- S1 <- diag(0, ncol(y))
      }
    } else if(model == "L"){
      if(hyper){
        if(is.null(prior$m1)){ m1 = rep(0, ncol(y)) } else { m1 = prior$m1 }
        if(is.null(prior$k1)){ k1 = 1 } else { k1 = prior$k1 }
        if(is.null(prior$lambda1)){ lambda1 = ncol(y) + 2 } else { lambda1 = prior$lambda1 }
        if(is.null(prior$Lambda1)){ Lambda1 = diag(1, ncol(y)) } else { Lambda1 = prior$Lambda1 }
        if(is.null(prior$Sigma0)){ Sigma0 = diag(1, ncol(y)) } else { Sigma0 = prior$Sigma0 }
        if(is.null(prior$n0)){ n0 = ncol(y) + 2 } else { n0 = prior$n0 }

        if(is.null(prior$S20)){ S20 = solve(rWishart(n = 1, Sigma = solve(Lambda1), df = lambda1)[,,1]) } else {S20 = prior$S20}
        if(is.null(prior$m0)){ m0 = as.vector(rnorm(ncol(y)) %*% t(chol(S20)) + m1)} else { m0 = prior$m0 }


      } else {
        if(is.null(prior$m0)){ m0 = rep(0, ncol(y)) } else { m0 = prior$m0 }
        if(is.null(prior$S20)){ S20 = diag(1, ncol(y)) } else { S20 = prior$S20 }
        if(is.null(prior$Sigma0)){ Sigma0 = diag(1, ncol(y)) } else { Sigma0 = prior$Sigma0 }
        if(is.null(prior$n0)){ n0 = ncol(y) + 2 } else { n0 = prior$n0 }

        m1 <- rep(0, ncol(y))
        k1 <- lambda1 <- 0
        Lambda1 <- diag(0, ncol(y))
      }
    } else if(model == "DLS"){

      if(hyper){
        if(is.null(prior$a0)){ a0 = rep(2, ncol(y)) } else { a0 = prior$a0 }
        if(is.null(prior$m1)){ m1 = rep(0, ncol(y)) } else { m1 = prior$m1 }
        if(is.null(prior$s21)){ s21 = rep(1, ncol(y)) } else { s21 = prior$s21 }
        if(is.null(prior$tau1)){ tau1 = rep(1, ncol(y)) } else { tau1 = prior$tau1 }
        if(is.null(prior$zeta1)){ zeta1 = rep(1, ncol(y)) } else { zeta1 = prior$zeta1 }
        if(is.null(prior$a1)){ a1 = rep(1, ncol(y)) } else { a1 = prior$a1 }
        if(is.null(prior$b1)){ b1 = rep(1, ncol(y)) } else { b1 = prior$b1 }

        if(is.null(prior$m0)){ m0 = apply(cbind(m1, s21), 1,
                                              function(x) rnorm(1, x[1], sqrt(x[2]))) } else { m0 = prior$m0 }
        if(is.null(prior$k0)){ k0 = apply(cbind(tau1,zeta1), 1, function(x) rgamma(1, x[1], x[2])) } else { k0 = prior$k0 }
        if(is.null(prior$b0)){ b0 = apply(cbind(a1,b1), 1, function(x) rgamma(1, x[1], x[2])) } else { b0 = prior$b0 }
      } else {
        if(is.null(prior$m0)){ m0 = rep(0, ncol(y)) } else { m0 = prior$m0 }
        if(is.null(prior$k0)){ k0 = rep(1, ncol(y)) } else { k0 = prior$k0 }
        if(is.null(prior$a0)){ a0 = rep(2, ncol(y)) } else { a0 = prior$a0 }
        if(is.null(prior$b0)){ b0 = rep(1, ncol(y)) } else { b0 = prior$b0 }

        m1 <- s21 <- tau1 <- zeta1 <- a1 <- b1 <- rep(0, ncol(y))
      }
    }

    # process parameters
    strength = ifelse(is.null(prior$strength), 1, prior$strength)
    discount = ifelse(is.null(prior$discount), 0, prior$discount)
    if(strength < - discount) stop("strength must be greater than -discount")

    # convert data to matrix
    y <- as.matrix(y)

    # estimate the model
    if(method == "ICS"){
      if(model == "LS"){
        est_model <- cICS_mv(y, grid_use, niter, nburn, m0, k0, Sigma0, n0, m1, S1, tau1, zeta1, n1, Sigma1,
                             strength, m_imp, nupd, out_param, mcmc_dens, discount, print_message, mean_dens, hyper)
      } else if(model == "L"){
        est_model <- cICS_mv_L(y, grid_use, niter, nburn, m0, S20, Sigma0, n0, m1, k1, lambda1, Lambda1,
                               strength, m_imp, nupd, out_param, mcmc_dens, discount, print_message, mean_dens, hyper)
      } else if(model == "DLS"){
        est_model <- cICS_mv_P(y, grid_use, niter, nburn, m0, k0, a0, b0, m1, s21, tau1, zeta1, a1, b1,
                               strength, m_imp, nupd, out_param, mcmc_dens, discount, print_message, mean_dens, hyper)
      }
    } else if(method == "SLI"){
      if(model == "LS"){
        est_model <- cSLI_mv(y, grid_use, niter, nburn, m0, k0, Sigma0, n0, m1, S1, tau1, zeta1, n1, Sigma1,
                             strength, nupd, out_param, mcmc_dens, discount, print_message, mean_dens, hyper)
      } else if(model == "L"){
        est_model <- cSLI_mv_L(y, grid_use, niter, nburn, m0, S20, Sigma0, n0, m1, k1, lambda1, Lambda1,
                               strength, nupd, out_param, mcmc_dens, discount, print_message, mean_dens, hyper)
      } else if(model == "DLS"){
        est_model <- cSLI_mv_P(y, grid_use, niter, nburn, m0, k0, a0, b0, m1, s21, tau1, zeta1, a1, b1,
                               strength, nupd, out_param, mcmc_dens, discount, print_message, mean_dens, hyper)
      }
    } else if(method == "MAR"){
      if(model == "LS"){
        est_model <- MAR_mv(y, grid_use, niter, nburn, m0, k0, Sigma0, n0, m1, S1, tau1, zeta1, n1, Sigma1,
                            strength, nupd, out_param, mcmc_dens, discount, print_message, mean_dens, hyper)
      } else if(model == "L"){
        est_model <- MAR_mv_L(y, grid_use, niter, nburn, m0, S20, Sigma0, n0, m1, k1, lambda1, Lambda1,
                              strength, nupd, out_param, mcmc_dens, discount, print_message, mean_dens, hyper)
      } else if(model == "DLS"){
        est_model <- MAR_mv_P(y, grid_use, niter, nburn, m0, k0, a0, b0, m1, s21, tau1, zeta1, a1, b1,
                              strength, nupd, out_param, mcmc_dens, discount, print_message, mean_dens, hyper)
      }
    }

    # return the results
    if(!isTRUE(out_param)){
      if(isTRUE(mcmc_dens)){
        result <- BNPdens(density = est_model$dens,
                          data = y,
                          grideval = grid_use,
                          clust = (est_model$clust + 1),
                          niter = niter,
                          nburn = nburn,
                          tot_time = est_model$time,
                          univariate = FALSE)
      }else{
        result <- BNPdens(data = y,
                          clust = (est_model$clust + 1),
                          niter = niter,
                          nburn = nburn,
                          tot_time = est_model$time,
                          univariate = FALSE)
      }
    } else {
      if(isTRUE(mcmc_dens)){
        result <- BNPdens(density = est_model$dens,
                          data = y,
                          grideval = grid_use,
                          clust = (est_model$clust + 1),
                          mean = est_model$mu,
                          sigma2 = est_model$s2,
                          probs = est_model$probs,
                          niter = niter,
                          nburn = nburn,
                          tot_time = est_model$time,
                          univariate = FALSE)
      }else{
        result <- BNPdens(data = y,
                          clust = (est_model$clust + 1),
                          mean = est_model$mu,
                          sigma2 = est_model$s2,
                          probs = est_model$probs,
                          niter = niter,
                          nburn = nburn,
                          tot_time = est_model$time,
                          univariate = FALSE)
      }
    }
  }
  return(result)
}
