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
#' }
#'
#' @param mcmc a list of MCMC arguments:
#' \itemize{
#'   \item \code{niter} (mandatory), number of iterations.
#'
#'   \item \code{nburn} (mandatory), number of iterations to discard as burn-in.
#'
#'   \item \code{method}, the MCMC sampling method to be used, options are \code{'ICS'}, \code{'MAR'} and \code{'SLI'} (default is \code{'ICS'}). See details.
#'
#'   \item \code{model}, the type of model to be fitted (default is \code{'LS'}). See details.
#'
#'   \item \code{nupd}, argument controlling the number of iterations to be displayed on screen: the function reports
#'   on standard output every time \code{nupd} new iterations have been carried out (default is \code{niter/10}).
#'
#'   \item \code{print_message}, control option. If equal to \code{TRUE}, the status is printed
#'   to standard output every \code{nupd} iterations (default is \code{TRUE}).
#'
#'   \item \code{m_imp}, number of generated values for the importance sampling step of \code{method = 'ICS'} (default is 10). See details.
#'
#'   \item \code{slice_type}, when \code{method = 'SLI'} it specifies the type of slice sampler. Options are \code{'DEP'} for dependent slice-efficient, and \code{'INDEP'} for independent slice-efficient (default is \code{'DEP'}). See details.
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
#'  mean and variance of the base measure on the location parameter (default are sample mean and sample variance of the data);
#'  \code{a0} and \code{b0} are shape and scale parameters of the inverse gamma prior on the common scale parameter
#'  (default are 2 and the sample variance of the data).
#'  If \code{hyper = TRUE},  optional hyperpriors on the base measure's parameters are added:
#'  specifically, \code{m1} and \code{k1} are the mean parameter and the scale factor defining the
#'  normal hyperprior on \code{m0} (default are the sample mean of the data and 1), and
#'  \code{a1} and \code{b1} are shape and rate parameters of the inverse gamma hyperprior on \code{b0}
#'  (default are 2 and the sample variance of the data). See details.
#'
#' \item If \code{model = 'LS'} (location-scale mixture) and \code{y} is univariate:
#'
#'   \code{m0} and \code{k0} are the mean parameter and the scale factor defining the normal base measure
#'   on the location parameter (default are the sample mean of the data and 1), and \code{a0} and \code{b0} are
#'   shape and scale parameters of the inverse gamma base measure on the scale parameter (default are 2 and the sample variance of the data).
#'  If \code{hyper = TRUE},  optional hyperpriors on the base measure's parameters are added:
#'  specifically, \code{m1} and \code{s21} are mean and variance parameters of the normal hyperprior on
#'  \code{m0} (default are the sample mean and the sample variance of the data);
#'  \code{tau1} and \code{zeta1} are shape and rate parameters of the gamma hyperprior on
#'  \code{k0} (default is 1 for both);
#'  \code{a1} and \code{b1} are shape and rate parameters of the gamma hyperprior on
#'  \code{b0}  (default are the sample variance of the data and 1). See details.
#'
#' \item If \code{model = 'L'} (location mixture) and \code{y} is multivariate (\code{p}-variate):
#'
#' \code{m0} and \code{S20} are
#'  mean and covariance of the base measure on the location parameter (default are the sample mean and the sample covariance of the data);
#'  \code{Sigma0} and \code{n0} are the parameters of the inverse Whishart prior on
#'   the common scale matrix (default are the sample covariance of the data and \code{p}+2).
#'  If \code{hyper = TRUE}, optional hyperpriors on the base measure's parameters are added:
#'  specifically, \code{m1} and \code{k1} are the mean parameter and the scale factor defining the
#'  normal hyperprior on \code{m0} (default are the sample mean of the data and 1), and
#'  \code{lambda} and \code{Lambda1} are the parameters (degrees of freedom and scale) of the inverse Wishart prior on \code{S20}
#'  (default are \code{p}+2 and the sample covariance of the data). See details.
#'
#' \item If \code{model = 'LS'} (location-scale mixture) and \code{y} is multivariate (\code{p}-variate):
#'
#' \code{m0} and \code{k0} are the mean parameter and the scale factor defining the normal base measure on the
#' location parameter (default are the sample mean of the data and 1), and
#' \code{n0} and \code{Sigma0} are the parameters (degrees of freedom and scale) of the inverse Wishart base measure on the location parameter
#' (default are \code{p}+2 and the sample covariance of the data).
#'  If \code{hyper = TRUE}, optional hyperpriors on the base measure's parameters are added:
#'  specifically, \code{m1} and \code{S1} are mean and covariance matrix parameters of the normal hyperprior on
#'  \code{m0} (default are the sample mean and the sample covariance of the data);
#'  \code{tau1} and \code{zeta1} are shape and rate parameters of the gamma hyperprior on
#'  \code{k0} (default is 1 for both);
#'  \code{n1} and \code{Sigma1} are the parameters (degrees of freedom and scale) of the Wishart prior for \code{Sigma0}
#'  (default are \code{p}+2 and the sample covariance of the data divided \code{p}+2). See details.
#'
#'
#' \item If \code{model = 'DLS'} (diagonal location-scale mixture):
#'
#' \code{m0} and \code{k0} are the mean vector parameter and the vector of scale factors defining the normal base measure
#' on the location parameter (default are the sample mean and a vector of ones),
#' and \code{a0} and \code{b0} are vectors of
#' shape and scale parameters defining the base measure on the scale parameters (default are a vector of twos and the diagonal
#' of the sample covariance of the data).
#' If \code{hyper = TRUE}, optional hyperpriors on the base measure's parameters are added:
#'  specifically, \code{m1} and \code{s21} are vectors of mean and variance parameters for the normal hyperpriors on the components of
#'  \code{m0} (default are the sample mean and the diagonal of the sample covariance of the data);
#'  \code{tau1} and \code{zeta1} are vectors of shape and rate parameters of the gamma hyperpriors on the components of
#'  \code{k0} (default is a vector of ones for both);
#'  \code{a1} and \code{b1} are vectors of shape and rate parameters of the gamma hyperpriors on the components of
#'  \code{b0}  (default is the diagonal of the sample covariance of the data and a vector of ones). See details.
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
#' The default method is the importance conditional sampler \code{'ICS'} (Canale et al. 2019). Other options are
#' the marginal sampler \code{'MAR'} (Neal, 2000) and the slice sampler \code{'SLI'} (Kalli et al. 2011).
#' The importance conditional sampler performs an importance sampling step when updating the values of
#' individual parameters \eqn{\theta}, which requires to sample \code{m_imp} values from a suitable
#' proposal. Large values of \code{m_imp} are known to improve the mixing of the chain
#' at the cost of increased running time (Canale et al. 2019). Two options are available for the slice sampler,
#' namely the dependent slice-efficient sampler (\code{slice_type = 'DEP'}), which is set as default, and the
#' independent slice-efficient sampler (\code{slice_type = 'INDEP'}) (Kalli et al. 2011). See Corradin et al. (to appear)
#' for more details.
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
#' arXiv preprint, arXiv:1906.08147
#'
#' Corradin, R., Canale, A., Nipoti, B. (2021), BNPmix: An R Package for Bayesian Nonparametric Modeling via Pitman-Yor Mixtures,
#' Journal of Statistical Software, 100, doi:10.18637/jss.v100.i15
#'
#' Kalli, M., Griffin, J. E., and Walker, S. G. (2011), Slice sampling mixture models.
#' Statistics and Computing 21, 93-105, doi:10.1007/s11222-009-9150-y
#'
#' Neal, R. M. (2000), Markov Chain Sampling Methods for Dirichlet Process Mixture Models,
#' Journal of Computational and Graphical Statistics 9, 249-265, doi:10.2307/1390653
#'

PYdensity <- function(y,
                       mcmc = list(),
                       prior = list(),
                       output = list()){

  if(!is.list(mcmc)) stop("mcmc must be a list")
  if(!is.list(prior)) stop("prior must be a list")
  if(!is.list(output)) stop("output must be a list")

  if(!is.null(mcmc$niter) && (!is.numeric(mcmc$niter) | (mcmc$niter<1))) stop("mcmc$niter must be a positive integer")
  if(!is.null(mcmc$nburn) && (!is.numeric(mcmc$nburn) | (mcmc$nburn<1)) & (mcmc$nburn>mcmc$niter)) stop("mcmc$nburn must be a positive integer less than niter")
  if(!is.null(mcmc$nupd) && (!is.numeric(mcmc$nupd)  | (mcmc$nupd<1))) stop("mcmc$nupd must be a positive integer")
  if(!is.null(mcmc$m_imp) && (!is.numeric(mcmc$m_imp) | (mcmc$m_imp<1))) stop("mcmc$m_imp must be a positive integer")
  if(!is.null(mcmc$print_message) & (!is.logical(mcmc$print_message))) stop("mcmc$print_message must be a logical value")
  if(!is.null(mcmc$hyper) & !is.logical(mcmc$hyper)) stop("mcmc$hyper must be a logical value")
  if(!is.null(output$grid) & !is.vector(output$grid) & !is.matrix(output$grid) & !is.data.frame(output$grid)) stop("wrong grid specification")

  if(is.vector(y)){
    p = 1
    if(!is.null(prior$m0) & !is.numeric(prior$m0)) stop("prior$m0 must be a numerical value")
    if(!is.null(prior$k0) && (!is.numeric(prior$k0) | (prior$k0<=0))) stop("prior$k0 must be a numerical positive value")
    if(!is.null(prior$a0) && (!is.numeric(prior$a0) | (prior$a0<=0))) stop("prior$a0 must be a numerical positive value")
    if(!is.null(prior$b0) && (!is.numeric(prior$b0) | (prior$b0<=0))) stop("prior$b0 must be a numerical positive value")
    if(!is.null(prior$s20) && (!is.numeric(prior$s20) | (prior$s20<=0))) stop("prior$s20 must be a numerical positive value")

    if(!is.null(prior$m1) & !is.numeric(prior$m1)) stop("prior$m1 must be a numerical value")
    if(!is.null(prior$s21) && (!is.numeric(prior$s21) | (prior$s21<=0))) stop("prior$s21 must be a numerical positive value")
    if(!is.null(prior$a1) && (!is.numeric(prior$a1) | (prior$a1<=0))) stop("prior$a1 must be a numerical positive value")
    if(!is.null(prior$b1) && (!is.numeric(prior$b1) | (prior$b1<=0))) stop("prior$b1 must be a numerical positive value")
    if(!is.null(prior$k1) && (!is.numeric(prior$k1) | (prior$k1<=0))) stop("prior$k1 must be a numerical positive value")
    if(!is.null(prior$tau1) && (!is.numeric(prior$tau1) | (prior$tau1<=0))) stop("prior$tau1 must be a numerical positive value")
    if(!is.null(prior$zeta1) && (!is.numeric(prior$zeta1) | (prior$zeta1<=0))) stop("prior$zeta1 must be a numerical positive value")

    # mandatory parameters
    if(is.null(mcmc$niter)) stop("Missing number of iterations")
    if(is.null(mcmc$nburn)) mcmc$nburn = 0

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
        grid_use = seq(from = min(y) - 0.1 * diff(range(y)), to = max(y) + 0.1 * diff(range(y)), length.out = 30)
      } else {
        if(length(dim(output$grid)) > 1) stop("Wrong grid dimension")
        grid_use <- as.vector(output$grid)
      }
    } else if (output$out_type == "CLUST"){
      mcmc_dens = FALSE
      grid_use = 0
    }

    if(!(model == "LS" | model == "L")) stop("Wrong model setting")
    slice_type <- mcmc$slice_type
    if(is.null(slice_type)){ slice_type <- "DEP"}
    if(!(slice_type == "DEP" | slice_type == "INDEP")) stop("Wrong mcmc$slice_type setting")
    if(!(method == "ICS" | method == "SLI" | method == "MAR")) stop("Wrong method setting")

    if(is.null(mcmc$wei_slice)){
      indep_sli = "DEFAULT"
    } else {
      indep_sli = "CUSTOM"
    }
    if(length(mcmc$wei_slice) > 2) stop("Wrong mcmc$wei_slice setting")

    hyper = ifelse(is.null(mcmc$hyper), TRUE, mcmc$hyper)
    # Check for different model parameters
    # if null, initialize default parameters
    if(model == "LS"){
      if(hyper){
        if(is.null(prior$a0)){ a0 = 2 } else { a0 = prior$a0 }
        if(is.null(prior$m1)){ m1 = mean(y) } else { m1 = prior$m1 }
        if(is.null(prior$s21)){ s21 = var(y) } else { s21 = prior$s21 }
        if(is.null(prior$tau1)){ tau1 = 1 } else { tau1 = prior$tau1 }
        if(is.null(prior$zeta1)){ zeta1 = 1 } else { zeta1 = prior$zeta1 }
        if(is.null(prior$a1)){ a1 = var(y) } else { a1 = prior$a1 }
        if(is.null(prior$b1)){ b1 = 1 } else { b1 = prior$b1 }

        if(is.null(prior$m0)){ m0 = rnorm(1, m1, sqrt(s21)) } else { m0 = prior$m0 }
        if(is.null(prior$k0)){ k0 = rgamma(1, tau1, zeta1) } else { k0 = prior$k0 }
        if(is.null(prior$b0)){ b0 = rgamma(1, a1, b1) } else { b0 = prior$b0 }
      } else {
        if(is.null(prior$m0)){ m0 = mean(y) } else { m0 = prior$m0 }
        if(is.null(prior$k0)){ k0 = 1 } else { k0 = prior$k0 }
        if(is.null(prior$a0)){ a0 = 2 } else { a0 = prior$a0 }
        if(is.null(prior$b0)){ b0 = var(y) } else { b0 = prior$b0 }

        m1 <- s21 <- tau1 <- zeta1 <- a1 <- b1 <- 0
      }
    } else if(model == "L"){
      if(hyper){
        if(is.null(prior$a0)){ a0 = 2 } else { a0 = prior$a0 }
        if(is.null(prior$b0)){ b0 = var(y) } else { b0 = prior$b0 }
        if(is.null(prior$m1)){ m1 = mean(y) } else { m1 = prior$m1 }
        if(is.null(prior$k1)){ k1 = 1 } else { k1 = prior$k1 }
        if(is.null(prior$a1)){ a1 = 2 } else { a1 = prior$a1 }
        if(is.null(prior$b1)){ b1 = var(y) } else { b1 = prior$b1 }

        if(is.null(prior$m0)){ m0 = rnorm(1, m1, sqrt(b1 / ((a1 - 1) * k1))) } else { m0 = prior$m0 }
        if(is.null(prior$s20)){ s20 = 1 / rgamma(1, a1, 1/b1) } else { s20 = prior$s20 }
      } else {
        if(is.null(prior$m0)){ m0 = mean(y) } else { m0 = prior$m0 }
        if(is.null(prior$s20)){ s20 = var(y) } else { s20 = prior$s20 }
        if(is.null(prior$a0)){ a0 = 2 } else { a0 = prior$a0 }
        if(is.null(prior$b0)){ b0 = var(y) } else { b0 = prior$b0 }

        m1 <- k1 <- a1 <- b1 <- 0
      }
    }

    # process parameters
    strength = ifelse(is.null(prior$strength), 1, prior$strength)
    discount = ifelse(is.null(prior$discount), 0, prior$discount)
    if(strength < - discount) stop("strength must be greater than -discount")
    if(is.null(mcmc$wei_slice)){
      mcmc$wei_slice <- c(strength, discount)
    }

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
    } else if(method == "SLI" & slice_type == "DEP"){

      # call the SLI univariate function
      if(model == "LS"){
        est_model <- cSLI(y, grid_use, niter, nburn, m0, k0, a0, b0, m1, s21, tau1, zeta1, a1, b1,
                          strength, mcmc$wei_slice[1], mcmc$wei_slice[2], nupd, out_param, mcmc_dens, discount, print_message, hyper, FALSE)
      } else if(model == "L"){
        est_model <- cSLI_L(y, grid_use, niter, nburn, m0, s20, a0, b0, m1, k1, a1, b1,
                            strength, mcmc$wei_slice[1], mcmc$wei_slice[2], nupd, out_param, mcmc_dens, discount, print_message, hyper, FALSE)
      }
    } else if(method == "SLI" & slice_type == "INDEP"){

      if(indep_sli == "DEFAULT"){
        # call the SLI univariate function
        if(model == "LS"){
          est_model <- cSLI(y, grid_use, niter, nburn, m0, k0, a0, b0, m1, s21, tau1, zeta1, a1, b1,
                            strength, mcmc$wei_slice[1], mcmc$wei_slice[2], nupd, out_param, mcmc_dens, discount, print_message, hyper, TRUE)
        } else if(model == "L"){
          est_model <- cSLI_L(y, grid_use, niter, nburn, m0, s20, a0, b0, m1, k1, a1, b1,
                              strength, mcmc$wei_slice[1], mcmc$wei_slice[2], nupd, out_param, mcmc_dens, discount, print_message, hyper, TRUE)
        }
      }

      if(indep_sli == "CUSTOM"){
        if(model == "LS"){
          est_model <- cSLI(y, grid_use, niter, nburn, m0, k0, a0, b0, m1, s21, tau1, zeta1, a1, b1,
                            strength, mcmc$wei_slice[1], mcmc$wei_slice[2], nupd, out_param, mcmc_dens, discount, print_message, hyper, TRUE)
        } else if(model == "L"){
          est_model <- cSLI_L(y, grid_use, niter, nburn, m0, s20, a0, b0, m1, k1, a1, b1,
                              strength, mcmc$wei_slice[1], mcmc$wei_slice[2], nupd, out_param, mcmc_dens, discount, print_message, hyper, TRUE)
        }
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
                          univariate = TRUE,
                          wvals = est_model$bound,
                          group_log = est_model$tdns)
      }else{
        result <- BNPdens(data = y,
                          clust = est_model$clust,
                          niter = niter,
                          nburn = nburn,
                          tot_time = est_model$time,
                          univariate = TRUE,
                          wvals = est_model$bound,
                          group_log = est_model$tdns)
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
                          univariate = TRUE,
                          wvals = est_model$bound,
                          group_log = est_model$tdns)
      }else{
        result <- BNPdens(data = y,
                          clust = est_model$clust,
                          mean = est_model$mu,
                          sigma2 = est_model$s2,
                          probs = est_model$probs,
                          niter = niter,
                          nburn = nburn,
                          tot_time = est_model$time,
                          univariate = TRUE,
                          wvals = est_model$bound,
                          group_log = est_model$tdns)
      }
    }
  } else if(!is.vector(y)){
    p = ncol(y)
    if(!is.null(prior$m0) && (!is.numeric(prior$m0) | length(prior$m0)!=p)) stop("prior$m0 must be a numerical vector of proper dimension")
    if(!is.null(prior$k0) && (!is.numeric(prior$k0))) stop("prior$k0 must be a numerical vector of proper dimension")
    if(!is.null(prior$a0) && (!is.numeric(prior$a0) | length(prior$a0)!=p)) stop("prior$a0 must be a numerical vector of proper dimension")
    if(!is.null(prior$b0) && (!is.numeric(prior$b0) | length(prior$b0)!=p)) stop("prior$b0 must be a numerical vector of proper dimension")
    if(!is.null(prior$a1) && (!is.numeric(prior$a1) | length(prior$a1)!=p)) stop("prior$a1 must be a numerical vector of proper dimension")
    if(!is.null(prior$b1) && (!is.numeric(prior$b1) | length(prior$b1)!=p)) stop("prior$b1 must be a numerical vector of proper dimension")
    if(!is.null(prior$k1) && (!is.numeric(prior$k1) | length(prior$k1)!=p)) stop("prior$k1 must be a numerical vector of proper dimension")
    if(!is.null(prior$Sigma0) && (!is.matrix(prior$Sigma0) | ncol(prior$Sigma0) != nrow(prior$Sigma0) | ncol(prior$Sigma0) !=p)) stop("prior$Sigma0 must be a square matrix of proper dimension")
    if(!is.null(prior$S20) && (!is.matrix(prior$S20) | ncol(prior$S0) != nrow(prior$S0) | ncol(prior$S0) !=p)) stop("prior$S20 must be a square matrix of proper dimension")
    if(!is.null(prior$n0) && (!is.numeric(prior$n0) | prior$n0<=(p+1) )) stop("prior$n0 must be a positive value grater than ncol(y) + 1")

    if(!is.null(prior$m1) && ( !is.vector(prior$m1) | length(prior$m1)!=p)) stop("prior$m1 must be a numerical vector of proper dimension")
    if(!is.null(prior$s21) && ( !is.vector(prior$s21) | length(prior$s21)!=p)) stop("prior$s21 must be a numerical vector of proper dimension")
    if(!is.null(prior$S1) && ( !is.matrix(prior$S1) )) stop("prior$S1 must be a square matrix of proper dimension")
    if(!is.null(prior$Sigma1) && ( !is.matrix(prior$Sigma1) | ncol(prior$Sigma1) != nrow(prior$Sigma1) | ncol(prior$Sigma1) !=p)) stop("prior$Sigma1 must be a square matrix of proper dimension")
    if(!is.null(prior$Lambda1) && ( !is.matrix(prior$Lambda1) | ncol(prior$Lambda1) != nrow(prior$Lambda1) | ncol(prior$Lambda1) !=p)) stop("prior$Lambda1 must be a square matrix of proper dimension")
    if(!is.null(prior$tau1) & !is.numeric(prior$tau1) & !is.vector(prior$tau1)) stop("prior$tau1 must be a numerical value or a vector")
    if(!is.null(prior$zeta1) & !is.numeric(prior$zeta1) & !is.vector(prior$zeta1)) stop("prior$zeta1 must be a numerical value or a vector")
    if(!is.null(prior$n1) && ( !is.numeric(prior$n1) | prior$n1<=(p-1) )) stop("prior$n1 must be a numerical value grater than ncol(y) - 1")
    if(!is.null(prior$lambda1) && ( !is.numeric(prior$lambda1))) stop("prior$lambda1 must be a numerical value")
    if(!is.null(prior$a1) && ( !is.vector(prior$a1) | length(prior$a1)!=p)) stop("prior$a1 must be a vector of proper dimension")
    if(!is.null(prior$b1) && ( !is.vector(prior$b1) | length(prior$b1)!=p)) stop("prior$b1 must be a vector of proper dimension")

    # Check the mandatory parameters
    if(is.null(mcmc$niter)) stop("Missing number of iterations")
    if(is.null(mcmc$nburn)) mcmc$nburn = 0

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
        if(p == 2){
          grid_use <- as.matrix(expand.grid(
            seq(from = min(y[,1]) - 0.1 * diff(range(y[,1])), to = max(y[,1]) + 0.1 * diff(range(y[,1])), length.out = 30),
            seq(from = min(y[,2]) - 0.1 * diff(range(y[,2])), to = max(y[,2]) + 0.1 * diff(range(y[,2])), length.out = 30)))
        } else if(p == 3){
          grid_use <- as.matrix(expand.grid(
            seq(from = min(y[,1]) - 0.1 * diff(range(y[,1])), to = max(y[,1]) + 0.1 * diff(range(y[,1])), length.out = 20),
            seq(from = min(y[,2]) - 0.1 * diff(range(y[,2])), to = max(y[,2]) + 0.1 * diff(range(y[,2])), length.out = 20),
            seq(from = min(y[,3]) - 0.1 * diff(range(y[,3])), to = max(y[,3]) + 0.1 * diff(range(y[,3])), length.out = 20)))
        } else {
          grid_use = matrix(0, nrow = 1, ncol = p)
        }
      } else {
        # check the grid
        if(ncol(output$grid) != ncol(y)) stop("The dimensions of grid and data not match")
        grid_use <- as.matrix(output$grid)
      }
    } else if (output$out_type == "MEAN"){
      mean_dens = TRUE
      mcmc_dens = TRUE
      if(is.null(output$grid)){
        if(p == 2){
          grid_use <- as.matrix(expand.grid(
            seq(from = min(y[,1]) - 0.1 * diff(range(y[,1])), to = max(y[,1]) + 0.1 * diff(range(y[,1])), length.out = 30),
            seq(from = min(y[,2]) - 0.1 * diff(range(y[,2])), to = max(y[,2]) + 0.1 * diff(range(y[,2])), length.out = 30)))
        } else if(p == 3){
          grid_use <- as.matrix(expand.grid(
            seq(from = min(y[,1]) - 0.1 * diff(range(y[,1])), to = max(y[,1]) + 0.1 * diff(range(y[,1])), length.out = 20),
            seq(from = min(y[,2]) - 0.1 * diff(range(y[,2])), to = max(y[,2]) + 0.1 * diff(range(y[,2])), length.out = 20),
            seq(from = min(y[,3]) - 0.1 * diff(range(y[,3])), to = max(y[,3]) + 0.1 * diff(range(y[,3])), length.out = 20)))
        } else {
          grid_use = matrix(0, nrow = 1, ncol = p)
        }
      } else {
        # check the grid
        if(ncol(output$grid) != ncol(y)) stop("The dimensions of grid and data not match")
        grid_use <- as.matrix(output$grid)
      }
    } else if (output$out_type == "CLUST"){
      mean_dens = FALSE
      mcmc_dens = FALSE
      grid_use = matrix(0, nrow = 1, ncol = p)
    }

    if(!(model == "LS" | model == "L" | model == "DLS")) stop("Wrong model setting")
    slice_type <- mcmc$slice_type
    if(is.null(slice_type)){ slice_type <- "DEP"}
    if(!(slice_type == "DEP" | slice_type == "INDEP")) stop("Wrong mcmc$slice_type setting")
    if(!(method == "ICS" | method == "SLI" | method == "MAR")) stop("Wrong method setting")

    if(is.null(mcmc$wei_slice)){
      indep_sli = "DEFAULT"
    } else {
      indep_sli = "CUSTOM"
    }
    if(length(mcmc$wei_slice) > 2) stop("Wrong mcmc$wei_slice setting")

    hyper = ifelse(is.null(mcmc$hyper), TRUE, mcmc$hyper)
    # Check for different model parameters
    if(model == "LS"){
      if(hyper){
        if(is.null(prior$m1)){ m1 = colMeans(y) } else { m1 = prior$m1 }
        if(is.null(prior$S1)){ S1 = var(y) } else { S1 = prior$S1 }
        if(is.null(prior$tau1)){ tau1 = 1 } else { tau1 = prior$tau1 }
        if(is.null(prior$zeta1)){ zeta1 = 1 } else { zeta1 = prior$zeta1 }
        if(is.null(prior$n1)){ n1 = ncol(y) + 2 } else { n1 = prior$n1 }
        if(is.null(prior$Sigma1)){ Sigma1 = var(y) } else { Sigma1 = prior$Sigma1 }
        if(is.null(prior$n0)){ n0 = ncol(y) + 2 } else { n0 = prior$n0 }

        if(is.null(prior$Sigma0)){ Sigma0 = rWishart(n = 1, Sigma = Sigma1, df = n1)[,,1] } else {Sigma0 = prior$Sigma0}
        if(is.null(prior$m0)){ m0 = as.vector(rnorm(ncol(y)) %*% t(chol(S1)) + m1)} else { m0 = prior$m0 }
        if(is.null(prior$k0)){ k0 = rgamma(1, tau1, 1 / zeta1)} else { k0 = prior$k0 }
      } else {
        if(is.null(prior$m0)){ m0 = colMeans(y) } else { m0 = prior$m0 }
        if(is.null(prior$k0)){ k0 = 1 } else { k0 = prior$k0 }
        if(is.null(prior$n0)){ n0 = ncol(y) + 2 } else { n0 = prior$n0 }
        if(is.null(prior$Sigma0)){ Sigma0 = var(y) / n0 } else { Sigma0 = prior$Sigma0 }

        m1 <- rep(0, ncol(y))
        n1 <- tau1 <- zeta1 <- 0
        Sigma1 <- S1 <- diag(0, ncol(y))
      }
    } else if(model == "L"){
      if(hyper){
        if(is.null(prior$m1)){ m1 = colMeans(y) } else { m1 = prior$m1 }
        if(is.null(prior$k1)){ k1 = 1 } else { k1 = prior$k1 }
        if(is.null(prior$lambda1)){ lambda1 = ncol(y) + 2 } else { lambda1 = prior$lambda1 }
        if(is.null(prior$Lambda1)){ Lambda1 = var(y) } else { Lambda1 = prior$Lambda1 }
        if(is.null(prior$n0)){ n0 = ncol(y) + 2 } else { n0 = prior$n0 }
        if(is.null(prior$Sigma0)){ Sigma0 = var(y) } else { Sigma0 = prior$Sigma0 }

        if(is.null(prior$S20)){ S20 = solve(rWishart(n = 1, Sigma = solve(Lambda1), df = lambda1)[,,1]) } else {S20 = prior$S20}
        if(is.null(prior$m0)){ m0 = as.vector(rnorm(ncol(y)) %*% t(chol(S20)) + m1)} else { m0 = prior$m0 }


      } else {
        if(is.null(prior$m0)){ m0 = colMeans(y) } else { m0 = prior$m0 }
        if(is.null(prior$S20)){ S20 = diag(1, ncol(y)) } else { S20 = prior$S20 }
        if(is.null(prior$Sigma0)){ Sigma0 = var(y) } else { Sigma0 = prior$Sigma0 }
        if(is.null(prior$n0)){ n0 = ncol(y) + 2 } else { n0 = prior$n0 }

        m1 <- rep(0, ncol(y))
        k1 <- lambda1 <- 0
        Lambda1 <- diag(0, ncol(y))
      }
    } else if(model == "DLS"){

      if(hyper){
        if(is.null(prior$a0)){ a0 = rep(2, ncol(y)) } else { a0 = prior$a0 }
        if(is.null(prior$m1)){ m1 = colMeans(y) } else { m1 = prior$m1 }
        if(is.null(prior$s21)){ s21 = diag(var(y)) } else { s21 = prior$s21 }
        if(is.null(prior$tau1)){ tau1 = rep(1, ncol(y)) } else { tau1 = prior$tau1 }
        if(is.null(prior$zeta1)){ zeta1 = rep(1, ncol(y)) } else { zeta1 = prior$zeta1 }
        if(is.null(prior$a1)){ a1 = diag(var(y)) } else { a1 = prior$a1 }
        if(is.null(prior$b1)){ b1 = rep(1, ncol(y)) } else { b1 = prior$b1 }

        if(is.null(prior$m0)){ m0 = apply(cbind(m1, s21), 1,
                                              function(x) rnorm(1, x[1], sqrt(x[2]))) } else { m0 = prior$m0 }
        if(is.null(prior$k0)){ k0 = apply(cbind(tau1,zeta1), 1, function(x) rgamma(1, x[1], x[2])) } else { k0 = prior$k0 }
        if(is.null(prior$b0)){ b0 = apply(cbind(a1,b1), 1, function(x) rgamma(1, x[1], x[2])) } else { b0 = prior$b0 }
      } else {
        if(is.null(prior$m0)){ m0 = colMeans(y) } else { m0 = prior$m0 }
        if(is.null(prior$k0)){ k0 = rep(1, ncol(y)) } else { k0 = prior$k0 }
        if(is.null(prior$a0)){ a0 = rep(2, ncol(y)) } else { a0 = prior$a0 }
        if(is.null(prior$b0)){ b0 = diag(var(y)) } else { b0 = prior$b0 }

        m1 <- s21 <- tau1 <- zeta1 <- a1 <- b1 <- rep(0, ncol(y))
      }
    }

    # process parameters
    strength = ifelse(is.null(prior$strength), 1, prior$strength)
    discount = ifelse(is.null(prior$discount), 0, prior$discount)
    if(strength < - discount) stop("strength must be greater than -discount")
    if(is.null(mcmc$wei_slice)){
      mcmc$wei_slice <- c(strength, discount)
    }

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
    } else if(method == "SLI" & slice_type == "DEP"){
      if(model == "LS"){
        est_model <- cSLI_mv(y, grid_use, niter, nburn, m0, k0, Sigma0, n0, m1, S1, tau1, zeta1, n1, Sigma1,
                             strength, mcmc$wei_slice[1], mcmc$wei_slice[2], nupd, out_param, mcmc_dens, discount, print_message, mean_dens, hyper, FALSE)
      } else if(model == "L"){
        est_model <- cSLI_mv_L(y, grid_use, niter, nburn, m0, S20, Sigma0, n0, m1, k1, lambda1, Lambda1,
                               strength, mcmc$wei_slice[1], mcmc$wei_slice[2], nupd, out_param, mcmc_dens, discount, print_message, mean_dens, hyper, FALSE)
      } else if(model == "DLS"){
        est_model <- cSLI_mv_P(y, grid_use, niter, nburn, m0, k0, a0, b0, m1, s21, tau1, zeta1, a1, b1,
                               strength, mcmc$wei_slice[1], mcmc$wei_slice[2], nupd, out_param, mcmc_dens, discount, print_message, mean_dens, hyper, FALSE)
      }
    } else if(method == "SLI" & slice_type == "INDEP"){
      if(indep_sli == "DEFAULT"){
        if(model == "LS"){
          est_model <- cSLI_mv(y, grid_use, niter, nburn, m0, k0, Sigma0, n0, m1, S1, tau1, zeta1, n1, Sigma1,
                               strength, mcmc$wei_slice[1], mcmc$wei_slice[2], nupd, out_param, mcmc_dens, discount, print_message, mean_dens, hyper, TRUE)
        } else if(model == "L"){
          est_model <- cSLI_mv_L(y, grid_use, niter, nburn, m0, S20, Sigma0, n0, m1, k1, lambda1, Lambda1,
                                 strength, mcmc$wei_slice[1], mcmc$wei_slice[2], nupd, out_param, mcmc_dens, discount, print_message, mean_dens, hyper, TRUE)
        } else if(model == "DLS"){
          est_model <- cSLI_mv_P(y, grid_use, niter, nburn, m0, k0, a0, b0, m1, s21, tau1, zeta1, a1, b1,
                                 strength, mcmc$wei_slice[1], mcmc$wei_slice[2], nupd, out_param, mcmc_dens, discount, print_message, mean_dens, hyper, TRUE)
        }
      }

      if(indep_sli == "CUSTOM"){
        if(model == "LS"){
          est_model <- cSLI_mv(y, grid_use, niter, nburn, m0, k0, Sigma0, n0, m1, S1, tau1, zeta1, n1, Sigma1,
                               strength, mcmc$wei_slice[1], mcmc$wei_slice[2], nupd, out_param, mcmc_dens,
                               discount, print_message, mean_dens, hyper, TRUE)
        } else if(model == "L"){
          est_model <- cSLI_mv_L(y, grid_use, niter, nburn, m0, S20, Sigma0, n0, m1, k1, lambda1, Lambda1,
                                 strength, mcmc$wei_slice[1], mcmc$wei_slice[2], nupd, out_param, mcmc_dens,
                                 discount, print_message, mean_dens, hyper, TRUE)
        } else if(model == "DLS"){
          est_model <- cSLI_mv_P(y, grid_use, niter, nburn, m0, k0, a0, b0, m1, s21, tau1, zeta1, a1, b1,
                                 strength, mcmc$wei_slice[1], mcmc$wei_slice[2], nupd, out_param, mcmc_dens,
                                 discount, print_message, mean_dens, hyper, TRUE)
        }
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
                          univariate = FALSE,
                          wvals = est_model$wvals)
      }else{
        result <- BNPdens(data = y,
                          clust = (est_model$clust + 1),
                          niter = niter,
                          nburn = nburn,
                          tot_time = est_model$time,
                          univariate = FALSE,
                          wvals = est_model$wvals)
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
                          univariate = FALSE,
                          wvals = est_model$wvals)
      }else{
        result <- BNPdens(data = y,
                          clust = (est_model$clust + 1),
                          mean = est_model$mu,
                          sigma2 = est_model$s2,
                          probs = est_model$probs,
                          niter = niter,
                          nburn = nburn,
                          tot_time = est_model$time,
                          univariate = FALSE,
                          wvals = est_model$wvals)
      }
    }
  }
  return(result)
}
