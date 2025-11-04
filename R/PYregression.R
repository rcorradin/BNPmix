#' @name PYregression
#' @export PYregression
#'
#' @title MCMC for Pitman-Yor mixture of Gaussian regressions
#'
#' @description The \code{PYregression} function generates a posterior sample
#' for mixtures of linear regression models inspired by the ANOVA-DDP model
#' introduced in De Iorio et al. (2004). See details below for model specification.
#'
#'@param y a vector of observations, univariate dependent variable;
#'@param x a matrix of observations, multivariate independent variable;

#'@param output list of posterior summaries:
#'
#'\itemize{
#
#'\item \code{grid_y}, a vector of points where to evaluate the estimated posterior mean density of
#'\code{y}, conditionally on each value of \code{x} in \code{grid_x};
#'
#'\item \code{grid_x}, a matrix of points where to evaluate the realization of the posterior conditional densities of
#'\code{y} given \code{x};
#' \item \code{out_type}, if \code{out_type = "FULL"}, the function returns the estimated partitions and the realizations of the posterior density for each iteration;
#' If \code{out_type = "MEAN"}, return the estimated partitions and the mean of the densities sampled at each iteration;
#' If \code{out_type = "CLUST"}, return the estimated partitions. Default \code{out_type = "FULL"};
#'\item \code{out_param}, if equal to \code{TRUE}, the function returns the draws of the kernel's
#'   parameters for each MCMC iteration, default is \code{FALSE}. See \code{value} for details.
#'
#'}
#'
#'@param mcmc a list of MCMC arguments:
#' \itemize{
#'   \item \code{niter} (mandatory), number of iterations.
#'
#'   \item \code{nburn} (mandatory), number of iterations to discard as burn-in.
#'
#'   \item \code{method}, the MCMC sampling method to be used. Options are \code{'ICS'}, \code{'MAR'} and \code{'SLI'} (default is \code{'ICS'}). See details.
#'
#'   \item \code{model} the type of model to be fitted (default is 'LS'). See details.
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
#'   \item \code{m_marginal}, number of generated values for the augmentation step needed, if \code{method = 'MAR'}, to implement Algorithm 8 of Neal, 2000. (Default is 100). See details.
#'
#'   \item \code{hyper}, if equal to \code{TRUE}, hyperprior distributions on the base measure's
#'   parameters are added, as specified in \code{prior} and explained in \code{details} (default is \code{TRUE}).
#' }
#'
#' @param prior a list giving the prior information. The list includes
#' \code{strength} and \code{discount}, the strength and discount parameters of the Pitman-Yor process
#' (default are \code{strength = 1} and \code{discount = 0}, the latter leading to the Dirichlet process).
#' The remaining parameters specify the base measure: \code{m0} and \code{S0} are
#'  the mean and covariance of normal base measure on the regression coefficients (default are a vector of zeroes, except for the first element equal
#'  to \code{mean(y)}, and a diagonal matrix with each element equal to 100);
#'  \code{a0} and \code{b0} are the shape and scale parameters of the inverse gamma base measure on the scale component
#'  (default are 2 and var(y)).
#'  If \code{hyper = TRUE},  optional hyperpriors on the base measure's parameters are added:
#'  specifically, \code{m1} and \code{k1} are the  mean parameter and scale factor defining the
#'  normal hyperprior on \code{m0} (default are a vector of zeroes, except for the first element equal
#'  to the sample mean of the dependent observed variable, and 1);
#'  \code{tau1} and \code{zeta1} are the shape and rate parameters of the gamma hyperprior on
#'  \code{b0} (default is 1 for both);
#'  \code{n1} and \code{S1} are the parameters (degrees of freedom and scale) of the Wishart prior for \code{S0}
#'  (default 4 and a diagonal matrix with each element equal to 100);  See details.
#'
#' @details
#' This function fits a Pitman-Yor process mixture of Gaussian linear regression models, i.e
#' \deqn{\tilde f(y) = \int \phi(y; x^T \beta, \sigma^2) \tilde p (d \beta, d \sigma^2)}{%
#'       \tilde f(y) = \int \phi(y; x^T \beta, \sigma^2) \tilde p (d \beta, d \sigma^2),}
#'       where \eqn{x} is a bivariate vector containing the dependent variable in \code{x} and a value of 1
#'        for the intercept term.
#' The mixing measure \eqn{\tilde p} has a Pitman-Yor process prior with strength \eqn{\vartheta},
#' discount parameter \eqn{\alpha}. The location model assume a base measures \eqn{P_0} specified as
#' \deqn{P_0(d \beta) = N(d \beta; m_0, S_0) .}{%
#'       P0(d \beta) = N(d \beta; m0, S0).}
#'  while the location-scale model assume a base measures \eqn{P_0} specified as
#' \deqn{P_0(d \beta, d \sigma^2) = N(d \beta; m_0, S_0) \times IGa(d \sigma^2; a_0, b_0).}{%
#'       P0(d \beta, d \sigma^2) = N(d \beta; m0, S0)  IG(d \sigma^2; a0, b0).}
#'  Optional hyperpriors complete the model specification:
#' \deqn{m_0 \sim N(m_1, S_0 / k_1 ),\quad S_0 \sim IW(\nu_1, S_1),\quad b_0 \sim G(\tau_1, \zeta_1).}{%
#' m_0 ~ N(m1, S0/k1),  S0 ~ IW(\nu1, S1),  b0 ~ G(\tau1, \zeta1).}
#'
#' \strong{Posterior simulation methods}
#'
#' This generic function implements three types of MCMC algorithms for posterior simulation.
#' The default method is the importance conditional sampler \code{'ICS'} (Canale et al. 2019). Other options are
#' the marginal sampler \code{'MAR'} (algorithm 8 of Neal, 2000) and the slice sampler \code{'SLI'} (Kalli et al. 2011).
#' The importance conditional sampler performs an importance sampling step when updating the values of
#' individual parameters \eqn{\theta}, which requires to sample \code{m_imp} values from a suitable
#' proposal. Large values of \code{m_imp} are known to improve the mixing of the posterior distribution
#' at the cost of increased running time (Canale et al. 2019). When updateing the individual parameter
#' \eqn{\theta}, Algorithm 8 of Neal, 2000, requires to sample \code{m_marginal} values from the base
#' measure. \code{m_marginal} can be chosen arbitrarily. Two options are available for the slice sampler,
#' namely the dependent slice-efficient sampler (\code{slice_type = 'DEP'}), which is set as default, and the
#' independent slice-efficient sampler (\code{slice_type = 'INDEP'}) (Kalli et al. 2011). See Corradin et al. (to appear)
#' for more details.
#'
#'
#' @return A \code{BNPdens} class object containing the estimated density and
#' the cluster allocations for each iterations. The output contains also the data and
#' the grids. If \code{out_param = TRUE} the output
#' contains also the kernel specific parameters for each iteration. If \code{mcmc_dens = TRUE}, the
#' function returns also a realization from the posterior density for each iteration.
#' If \code{mean_dens = TRUE}, the output contains just the mean of the densities sampled at each iteration.
#' The output retuns also the number of iterations,
#' the number of burn-in iterations, the computational time and the type of model.
#'
#'
#' @references
#'
#' Canale, A., Corradin, R., Nipoti, B. (2019), Importance conditional sampling for Bayesian nonparametric mixtures,
#' arXiv preprint, arXiv:1906.08147
#'
#' Corradin, R., Canale, A., Nipoti, B. (2021), BNPmix: An R Package for Bayesian Nonparametric Modeling via Pitman-Yor Mixtures,
#' Journal of Statistical Software, doi:10.18637/jss.v100.i15
#'
#' De Iorio, M., Mueller, P., Rosner, G.L., and MacEachern, S. (2004), An ANOVA Model for Dependent Random Measures,
#' Journal of the American Statistical Association 99, 205-215, doi:10.1198/016214504000000205
#'
#' Kalli, M., Griffin, J. E., and Walker, S. G. (2011), Slice sampling mixture models.
#' Statistics and Computing 21, 93-105, doi:10.1007/s11222-009-9150-y
#'
#' Neal, R. M. (2000), Markov Chain Sampling Methods for Dirichlet Process Mixture Models,
#' Journal of Computational and Graphical Statistics 9, 249-265, doi:10.2307/1390653
#'
#'
#' @examples
#' x_toy <- c(rnorm(100, 3, 1), rnorm(100, 3, 1))
#' y_toy <- c(x_toy[1:100] * 2 + 1, x_toy[101:200] * 6 + 1) + rnorm(200, 0, 1)
#' grid_x <- c(0, 1, 2, 3, 4, 5)
#' grid_y <- seq(0, 35, length.out = 50)
#' est_model <- PYregression(y = y_toy, x = x_toy,
#' mcmc = list(niter = 200, nburn = 100),
#' output = list(grid_x = grid_x, grid_y = grid_y))
#' summary(est_model)
#' plot(est_model)
#'

PYregression <- function(y, x,
                          mcmc = list(),
                          prior = list(),
                          output = list()){

  if(!is.vector(y)) stop("Wrong dimensions: y need to be a vector")
  if(!(is.data.frame(x) | is.matrix(x) | is.vector(x))) stop("Wrong dimensions: x need to be a vector or a matrix")

  if(!is.list(mcmc)) stop("mcmc must be a list")
  if(!is.list(prior)) stop("prior must be a list")
  if(!is.list(output)) stop("output must be a list")

  if(!is.null(mcmc$niter) && (!is.numeric(mcmc$niter) | (mcmc$niter<1))) stop("mcmc$niter must be a positive integer")
  if(!is.null(mcmc$nburn) && (!is.numeric(mcmc$nburn) | (mcmc$nburn<1)) & (mcmc$nburn>mcmc$niter)) stop("mcmc$nburn must be a positive integer less than niter")
  if(!is.null(mcmc$nupd) && (!is.numeric(mcmc$nupd)  | (mcmc$nupd<1))) stop("mcmc$nupd must be a positive integer")
  if(!is.null(mcmc$m_imp) && (!is.numeric(mcmc$m_imp) | (mcmc$m_imp<1))) stop("mcmc$m_imp must be a positive integer")
  if(!is.null(mcmc$print_message) & (!is.logical(mcmc$print_message))) stop("mcmc$print_message must be a logical value")
  if(!is.null(mcmc$hyper) & !is.logical(mcmc$hyper)) stop("mcmc$hyper must be a logical value")
  if(!is.null(mcmc$m_marginal) & !is.numeric(mcmc$m_marginal)) stop("mcmc$m_marginal must be a numerical value")

  if(!is.null(prior$m0) & !is.vector(prior$m0)) stop("prior$m0 must be a vector of size 2")
  if(!is.null(prior$S0) && (!is.matrix(prior$S0) | ncol(prior$S0) != nrow(prior$S0) | ncol(prior$S0) !=2) ) stop("prior$S0 must be a square matrix of dimension 2")
  if(!is.null(prior$a0) && (!is.numeric(prior$a0) | (prior$a0<1))) stop("prior$a0 must be a positive value")
  if(!is.null(prior$b0) && (!is.numeric(prior$b0) | (prior$b0<1))) stop("prior$n0 must be a positive value")
  if(!is.null(prior$m1) & !is.vector(prior$m1)) stop("prior$m1 must be a vector")
  if(!is.null(prior$S1) & !is.matrix(prior$S1)) stop("prior$S1 must be a matrix")
  if(!is.null(prior$k1) & !is.numeric(prior$k1)) stop("prior$k1 must be a numerical value")
  if(!is.null(prior$tau1) && (!is.numeric(prior$tau1) | (prior$tau1<1)))  stop("prior$tau1 must be a positive value")
  if(!is.null(prior$zeta1) && (!is.numeric(prior$zeta1) | (prior$zeta1<1)))  stop("prior$zeta1 must be a positive value")
  if(!is.null(prior$n1) & !is.numeric(prior$n1)) stop("prior$n1 must be a numerical value")
  if(!is.null(prior$strength) & !is.numeric(prior$strength)) stop("prior$strength must be a numerical value")
  if(!is.null(prior$discount) & !is.numeric(prior$discount)) stop("prior$discount must be a numerical value")

  # mandatory parameters
  if(is.null(mcmc$niter)) stop("Missing number of iterations")
  if(is.null(mcmc$nburn)) mcmc$nburn = 0

  # add variable for the intercept
  if(length(dim(x)) <= 1){
    x <- as.matrix(cbind(rep(1, length(x)), x))
  } else {
    x <- as.matrix(cbind(rep(1, nrow(x)), x))
  }
  d <- ncol(x)

  # if mcmc misses some parts, add default
  niter = mcmc$niter
  nburn = mcmc$nburn
  method = ifelse(is.null(mcmc$method), "ICS", mcmc$method)
  model = ifelse(is.null(mcmc$model), "LS", mcmc$model)
  nupd = ifelse(is.null(mcmc$nupd), round(niter / 10), mcmc$nupd)
  print_message = ifelse(is.null(mcmc$print_message), TRUE, mcmc$print_message)
  m_imp = ifelse(is.null(mcmc$m_imp), 10, mcmc$m_imp)
  m_marginal = ifelse(is.null(mcmc$m_marginal), 100, mcmc$m_marginal)

  # output
  out_param = ifelse(is.null(output$out_param), FALSE, output$out_param)
  output$out = ifelse(is.null(output$out), "FULL", output$out)
  if(output$out == "FULL"){
    mean_dens = FALSE
    mcmc_dens = TRUE
    if(is.null(output$grid_y) & is.null(output$grid_x)){
      grid_y = seq(from = min(y) - 0.1 * diff(range(y)), to = max(y) + 0.1 * diff(range(y)), length.out = 30)
      grid_x = as.matrix(cbind(rep(1,4), apply(x, 2, function(z) seq(from = min(z) - 0.1 * diff(range(z)), to = max(z) + 0.1 * diff(range(z)), length.out = 4))[,-1]))
    } else if(is.null(output$grid_y) & !is.null(output$grid_x)){
      grid_y = seq(from = min(y) - 0.1 * diff(range(y)), to = max(y) + 0.1 * diff(range(y)), length.out = 30)
      if(ncol(x) == 2){
        grid_x = as.matrix(cbind(rep(1, length(output$grid_x)), output$grid_x))
      } else{
        grid_x = as.matrix(cbind(rep(1, nrow(output$grid_x)), output$grid_x))
      }
    } else if(!is.null(output$grid_y) & is.null(output$grid_x)){
      grid_y = output$grid_y
      grid_x = as.matrix(cbind(rep(1,4), apply(x, 2, function(z) seq(from = min(z) - 0.1 * diff(range(z)), to = max(z) + 0.1 * diff(range(z)), length.out = 4))[,-1]))
    } else {
      if(ncol(x) == 2){
        grid_x = as.matrix(cbind(rep(1, length(output$grid_x)), output$grid_x))
      } else{
        grid_x = as.matrix(cbind(rep(1, nrow(output$grid_x)), output$grid_x))
      }
      grid_y = output$grid_y
    }
  } else if (output$out == "MEAN"){
    mean_dens = TRUE
    mcmc_dens = TRUE
    if(is.null(output$grid_y) & is.null(output$grid_x)){
      grid_y = seq(from = min(y) - 0.1 * diff(range(y)), to = max(y) + 0.1 * diff(range(y)), length.out = 30)
      grid_x = as.matrix(cbind(rep(1,4), apply(x, 2, function(z) seq(from = min(z) - 0.1 * diff(range(z)), to = max(z) + 0.1 * diff(range(z)), length.out = 4))[,-1]))
    } else if(is.null(output$grid_y) & !is.null(output$grid_x)){
      grid_y = seq(from = min(y) - 0.1 * diff(range(y)), to = max(y) + 0.1 * diff(range(y)), length.out = 30)
      if(ncol(x) == 2){
        grid_x = as.matrix(cbind(rep(1, length(output$grid_x)), output$grid_x))
      } else{
        grid_x = as.matrix(cbind(rep(1, nrow(output$grid_x)), output$grid_x))
      }
    } else if(!is.null(output$grid_y) & is.null(output$grid_x)){
      grid_y = output$grid_y
      grid_x = as.matrix(cbind(rep(1,4), apply(x, 2, function(z) seq(from = min(z) - 0.1 * diff(range(z)), to = max(z) + 0.1 * diff(range(z)), length.out = 4))[,-1]))
    } else {
      if(ncol(x) == 2){
        grid_x = as.matrix(cbind(rep(1, length(output$grid_x)), output$grid_x))
      } else{
        grid_x = as.matrix(cbind(rep(1, nrow(output$grid_x)), output$grid_x))
      }
      grid_y = output$grid_y
    }
  } else if (output$out == "CLUST"){
    mean_dens = FALSE
    mcmc_dens = FALSE
    grid_y = seq(from = min(y) - 0.1 * diff(range(y)), to = max(y) + 0.1 * diff(range(y)), length.out = 30)
    grid_x = as.matrix(cbind(rep(1,4), apply(x, 2, function(z) seq(from = min(z) - 0.1 * diff(range(z)), to = max(z) + 0.1 * diff(range(z)), length.out = 4))[,-1]))
  }

  # mcmc_dens = ifelse(is.null(output$mcmc_dens), TRUE, output$mcmc_dens)
  # mean_dens = ifelse(is.null(output$mean_dens), FALSE, output$mean_dens)
  # if(is.null(output$grid_y)){
  #   mcmc_dens = FALSE
  #   grid_y = 0
  #   grid_x = matrix(c(0,0), ncol = 2)
  # } else {
  #   grid_x <- cbind(rep(1, length(output$grid_x)), output$grid_x)
  #   grid_y = output$grid_y
  # }

  slice_type <- mcmc$slice_type
  if(is.null(slice_type)){ slice_type <- "DEP"}
  if(!(slice_type == "DEP" | slice_type == "INDEP")) stop("Wrong mcmc$slice_type setting")
  if(!(method == "ICS" | method == "SLI" | method == "MAR")) stop("Wrong method setting")
  if(!(model == "LS" | model == "L")) stop("Wrong model setting")
  hyper = ifelse(is.null(mcmc$hyper), TRUE, mcmc$hyper)

  if(is.null(mcmc$wei_slice)){
    indep_sli = "DEFAULT"
  } else {
    indep_sli = "CUSTOM"
  }
  if(length(mcmc$wei_slice) > 2) stop("Wrong mcmc$wei_slice setting")

  # if null, initialize default parameters
  if(hyper){
    if(is.null(prior$a0)){ a0 = 2 } else { a0 = prior$a0 }
    if(is.null(prior$m1)){ m1 = c(mean(y), rep(0, d - 1)) } else { m1 = prior$m1 }
    if(is.null(prior$k1)){ k1 = 1 } else { k1 = prior$k1 }
    if(is.null(prior$tau1)){ tau1 = var(y) } else { tau1 = prior$tau1 }
    if(is.null(prior$zeta1)){ zeta1 = 1 } else { zeta1 = prior$zeta1 }
    if(is.null(prior$n1)){ n1 = d + 2 } else { n1 = prior$n1 }
    if(is.null(prior$S1)){ S1 = diag(100, d) } else { S1 = prior$S1 }

    if(is.null(prior$S0)){ S0 = solve(rWishart(n = 1, Sigma = solve(S1), df = n1)[,,1]) } else { S0 = prior$S0 }
    if(is.null(prior$m0)){ m0 = as.vector(rnorm(d) %*% (t(chol(S0)) / k1) + m1) } else { m0 = prior$m0 }
    if(is.null(prior$b0)){ b0 = rgamma(1, tau1, zeta1) } else { b0 = prior$b0 }
  } else {
    if(is.null(prior$m0)){ m0 = c(mean(y), rep(0, d - 1))  } else { m0 = prior$m0 }
    if(is.null(prior$S0)){ S0 = diag(100, d) } else { S0 = prior$S0 }
    if(is.null(prior$a0)){ a0 = 2 } else { a0 = prior$a0 }
    if(is.null(prior$b0)){ b0 = var(y) } else { b0 = prior$b0 }

    m1 <- rep(0, d)
    k1 <- n1 <- 1
    tau1 <- zeta1 <- 0
    S1 <- diag(1, d)
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
    if(model == "LS"){
      est_model <- cICS_mv_MKR(y, x, grid_y, grid_x, niter, nburn, m0, S0, a0, b0,
                               m1, k1, n1, S1, tau1, zeta1, strength, m_imp, nupd, out_param,
                               mcmc_dens, discount, print_message, mean_dens, hyper)
    } else if(model == "L"){
      est_model <- cICS_mv_MKR_L(y, x, grid_y, grid_x, niter, nburn, m0, S0, a0, b0,
                                 m1, k1, n1, S1, strength, m_imp, nupd, out_param,
                                 mcmc_dens, discount, print_message, mean_dens, hyper)
    }
  } else if(method == "SLI" & slice_type == "INDEP"){
    if(model == "LS"){
      est_model <- cSLI_mv_MKR(y, x, grid_y, grid_x, niter, nburn, m0, S0, a0, b0,
                               m1, k1, n1, S1, tau1, zeta1, strength, mcmc$wei_slice[1], mcmc$wei_slice[2], nupd, out_param,
                               mcmc_dens, discount, print_message, mean_dens, hyper, TRUE)
    } else if(model == "L"){
      est_model <- cSLI_mv_MKR_L(y, x, grid_y, grid_x, niter, nburn, m0, S0, a0, b0,
                                 m1, k1, n1, S1, strength, mcmc$wei_slice[1], mcmc$wei_slice[2], nupd, out_param,
                                 mcmc_dens, discount, print_message, mean_dens, hyper, TRUE)
    }
  } else if(method == "SLI" & slice_type == "DEP"){
    if(model == "LS"){
      est_model <- cSLI_mv_MKR(y, x, grid_y, grid_x, niter, nburn, m0, S0, a0, b0,
                               m1, k1, n1, S1, tau1, zeta1, strength, mcmc$wei_slice[1], mcmc$wei_slice[2], nupd, out_param,
                               mcmc_dens, discount, print_message, mean_dens, hyper, FALSE)
    } else if(model == "L"){
      est_model <- cSLI_mv_MKR_L(y, x, grid_y, grid_x, niter, nburn, m0, S0, a0, b0,
                                 m1, k1, n1, S1, strength, mcmc$wei_slice[1], mcmc$wei_slice[2], nupd, out_param,
                                 mcmc_dens, discount, print_message, mean_dens, hyper, FALSE)
    }
  } else if(indep_sli == "CUSTOM"){
    if(model == "LS"){
      est_model <- cSLI_mv_MKR(y, x, grid_y, grid_x, niter, nburn, m0, S0, a0, b0,
                               m1, k1, n1, S1, tau1, zeta1, strength, mcmc$wei_slice[1], mcmc$wei_slice[2], nupd, out_param,
                               mcmc_dens, discount, print_message, mean_dens, hyper, FALSE)
    } else if(model == "L"){
      est_model <- cSLI_mv_MKR_L(y, x, grid_y, grid_x, niter, nburn, m0, S0, a0, b0,
                                 m1, k1, n1, S1, strength, mcmc$wei_slice[1], mcmc$wei_slice[2], nupd, out_param,
                                 mcmc_dens, discount, print_message, mean_dens, hyper, FALSE)
    }
  } else if(method == "MAR"){
    if(model == "LS"){
      est_model <- MAR_mv_MKR(y, x, grid_y, grid_x, niter, nburn, m0, S0, a0, b0,
                              m1, k1, n1, S1, tau1, zeta1, strength, m_marginal, nupd, out_param,
                              mcmc_dens, discount, print_message, mean_dens, hyper)
    } else if(model == "L"){
      est_model <- MAR_mv_MKR_L(y, x, grid_y, grid_x, niter, nburn, m0, S0, a0, b0,
                                m1, k1, n1, S1, strength, m_marginal, nupd, out_param,
                                mcmc_dens, discount, print_message, mean_dens, hyper)
    }
  }

  # return the results
  if(!isTRUE(out_param)){
    if(isTRUE(mcmc_dens)){
      result <- BNPdens(density = est_model$dens,
                        data = cbind(y, x[,-1]),
                        grid_y = grid_y,
                        grid_x = grid_x,
                        clust = est_model$clust,
                        niter = niter,
                        nburn = nburn,
                        tot_time = est_model$time,
                        regression = TRUE)
    }else{
      result <- BNPdens(clust = est_model$clust,
                        data = cbind(y, x[,-1]),
                        niter = niter,
                        nburn = nburn,
                        tot_time = est_model$time,
                        regression = TRUE)
    }
  } else {
    if(isTRUE(mcmc_dens)){
      result <- BNPdens(density = est_model$dens,
                        data = cbind(y, x[,-1]),
                        grid_y = grid_y,
                        grid_x = grid_x,
                        clust = est_model$clust,
                        beta = est_model$beta,
                        sigma2 = est_model$sigma2,
                        probs = est_model$probs,
                        niter = niter,
                        nburn = nburn,
                        tot_time = est_model$time,
                        regression = TRUE)
    }else{
      result <- BNPdens(clust = est_model$clust,
                        data = cbind(y, x[,-1]),
                        beta = est_model$beta,
                        sigma2 = est_model$sigma2,
                        probs = est_model$probs,
                        niter = niter,
                        nburn = nburn,
                        tot_time = est_model$time,
                        regression = TRUE)
    }
  }
  return(result)
}
