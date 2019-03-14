#' @name condMCMC
#' @export condMCMC
#'
#' @title MCMC for univariate Pitman-Yor mixtures
#' @description The condMCMC function estimate a univariate Pitman-Yor process mixture model
#' with Gaussian kernel. Three possible sampling strategies: importance conditional sampler,
#' slice sampler and marginal sampler.
#'
#' The models are of the form \deqn{\tilde f(x) = \int k(x, \theta) \tilde p (d \theta)} where
#' \eqn{k(x, \theta)} is an univariate gaussian kernel function, \eqn{\tilde p} is distributed as a Pitman-Yor
#' process with total mass \eqn{\vartheta}, discount parameter \eqn{\sigma} and normal-inverse-gamma base measure \eqn{P_0}, i.e.
#' \deqn{P_0 \sim N(\mu; m_0, k_0 \sigma^2) \times IG(\sigma^2; a_0, b_0).}
#'
#' @param data A dataset (vector).
#' @param grid A grid to evaluate the estimated density (vector).
#' @param niter Number of iterations to estimate the model.
#' @param nburn Number of burn-in iterations.
#' @param m0 Mean of the distribution of location component of the base measure.
#' @param k0 Tuning parameter of the location component variance of the base measure.
#' @param a0 Shap parameter of the scale component distribution of the base measure.
#' @param b0 Rate parameter of the scale component distribution of the base measure.
#' @param mass Mass parameter of the Dirichlet process.
#' @param method Different methods to estimate the model. Possible vaule "ICS", "SLI" or "MAR",
#' importance conditional sampler/slice sampler/marginal sampler, default "ICS".
#' @param napprox Number of approximating value for the conditioning distribution via "ICS" method, default 100.
#' @param nupd How frequently show the curren state of the estimation (number of iterations) - default 1000.
#' @param out_param If TRUE, save the parameters for each iteration, default FALSE.
#' @param out_dens If TRUE, return also the estimated density, default TRUE.
#' @param process Dirichlet process ("DP") or Pitman-Yor process ("PY"), default "DP".
#' @param sigma_PY Discount parameter of the Pitman-Yor process, default 0.
#' @param print_message If TRUE print the status of the estimation, default TRUE.
#'
#' @return A modCond class object contain the estimated density for each iterations,
#' the allocations for each iterations. If out_param is TRUE, also the parameters for each iteration.
#'
#' @examples
#' data_toy <- c(rnorm(100, -3, 1), rnorm(100, 3, 1))
#' grid <- seq(-7, 7, length.out = 50)
#' est_model <- condMCMC(data = data_toy, grid = grid, niter = 1000,
#'                       nburn = 100, napprox = 100, nupd = 100)
#' summary(est_model)
#' plot(est_model)
#'

condMCMC <- function(data, grid = NULL, niter, nburn, m0 = NULL, k0 = NULL,
                     a0 = NULL, b0 = NULL, mass = 1, method = "ICS",
                     napprox = 100, nupd = 1000, out_param = F, out_dens= TRUE,
                     process = "DP", sigma_PY = 0, print_message = TRUE){

  if(is.null(grid)) grid <- 0
  if(process == "DP") sigma_PY <- 0

  grid_use <- as.vector(grid)
  if(length(dim(data)) > 1){
    stop("The dataset must be a vector")
  }
  if(is.null(m0)) m0 <- 0
  if(is.null(k0)) k0 <- 1
  if(is.null(a0)) a0 <- 2
  if(is.null(b0)) b0 <- 1

  if(method == "ICS"){

    # call the ICS univariate function
    est_model <- cICS(data, grid_use, niter, nburn, m0, k0, a0, b0, mass,
                      napprox, nupd, out_param, out_dens, sigma_PY, print_message)
  } else if(method == "SLI"){

    # call the SLI univariate function
    est_model <- cSLI(data, grid_use, niter, nburn, m0, k0, a0, b0,
                      mass, nupd, out_param, out_dens, sigma_PY, print_message)
  } else if(method == "MAR"){

    # call the MAR univariate function
    est_model <- MAR(data, grid_use, niter, nburn, m0, k0, a0, b0,
                     mass, nupd, out_param, out_dens, sigma_PY, print_message)
  }

  if(!isTRUE(out_param)){
    if(isTRUE(out_dens)){
      output <- modCond(density = est_model$dens,
                        grideval = grid_use,
                        clust = est_model$clust,
                        niter = niter,
                        nburn = nburn,
                        nnew = as.vector(est_model$newval),
                        tot_time = est_model$time,
                        univariate = TRUE)
    }else{
      output <- modCond(clust = est_model$clust,
                        niter = niter,
                        nburn = nburn,
                        nnew = as.vector(est_model$newval),
                        tot_time = est_model$time,
                        univariate = TRUE)
    }
  } else {
    if(isTRUE(out_dens)){
      output <- modCond(density = est_model$dens,
                        grideval = grid_use,
                        clust = est_model$clust,
                        mean = est_model$mu,
                        sigma2 = est_model$s2,
                        probs = est_model$probs,
                        niter = niter,
                        nburn = nburn,
                        nnew = as.vector(est_model$newval),
                        tot_time = est_model$time,
                        univariate = TRUE)
    }else{
      output <- modCond(clust = est_model$clust,
                        mean = est_model$mu,
                        sigma2 = est_model$s2,
                        probs = est_model$probs,
                        niter = niter,
                        nburn = nburn,
                        nnew = as.vector(est_model$newval),
                        tot_time = est_model$time,
                        univariate = TRUE)
    }
  }


  return(output)
}
