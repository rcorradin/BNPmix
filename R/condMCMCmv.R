#' @name condMCMCmv
#' @export condMCMCmv
#'
#' @title MCMC for multivariate Pitman-Yor mixtures
#' @description The condMCMCmv function estimate a multivariate Pitman-Yor process mixture model with
#' Gaussian kernel. Three possible sampling strategies: importance conditional sampler,
#' slice sampler and marginal sampler.
#'
#' The models are of the form \deqn{\tilde f(\mathbf x) = \int k(\mathbf x, \mathbf \theta) \tilde p (d \mathbf \theta)} where
#' \eqn{k(\mathbf x, \mathbf \theta)} is a multivariate gaussian kernel function, \eqn{\tilde p} is distributed as a Pitman-Yor
#' process with total mass \eqn{\vartheta}, discount parameter \eqn{\sigma} and normal-inverse-wishart base measure \eqn{P_0}, i.e.
#' \deqn{P_0 \sim N(\mathbf \mu; \mathbf m_0, k_0 \Sigma) \times IW(\Sigma; n_0, S_0).}
#'
#' @param data A dataset (matrix).
#' @param grid A grid to evaluate the estimated density (matrix).
#' @param niter Number of iterations to estimate the model.
#' @param nburn Number of burn-in iterations.
#' @param m0 Mean of the distribution of location component of the base measure.
#' @param k0 Tuning parameter of the location component variance of the base measure.
#' @param S0 Matrix of the Inverse Wishart distribution of the scale component.
#' @param n0 Degree of freedom of the Inverse Wishart distribution of the scale component.
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
#' @param light_dens If TRUE return only the mean of the estimated densities, default TRUE.
#'
#' @return A modCondMv class object contain the estimated density for each iterations,
#' the allocations for each iterations. If out_param is TRUE, also the parameters for each iteration.
#'
#' @examples
#' data_toy <- cbind(c(rnorm(100, -3, 1), rnorm(100, 3, 1)),
#'                   c(rnorm(100, -3, 1), rnorm(100, 3, 1)))
#' grid <- expand.grid(seq(-7, 7, length.out = 50),
#'                     seq(-7, 7, length.out = 50))
#' est_model <- condMCMCmv(data = data_toy, grid = grid, niter = 1000,
#'                        nburn = 100, napprox = 100, nupd = 100)
#' summary(est_model)
#' plot(est_model)
#'

condMCMCmv <- function(data, grid = NULL, niter, nburn,  m0 = NULL, k0 = NULL,
                       S0 = NULL, n0 = NULL, mass = 1, method = "ICS",
                       napprox = 100, nupd = 1000, out_param = F, out_dens = TRUE,
                       process = "DP", sigma_PY = 0, print_message = TRUE,
                       light_dens = TRUE){

  if(is.null(grid)) grid <- matrix(0, ncol = ncol(data))
  if(process == "DP") sigma_PY <- 0

  grid_use <- as.matrix(grid)
  if(length(dim(data)) == 1){
    stop("Wrong call, use the univariate one")
  }

  if(ncol(grid) != ncol(data)){
    stop("The dimensions of grid and data are not matching")
  }

  if(is.null(m0)) m0 <- rep(0, ncol(data))
  if(is.null(k0)) k0 <- 1
  if(is.null(S0)) S0 <- diag(1, ncol(data))
  if(is.null(n0)) n0 <- ncol(data) + 2

  if(method == "ICS"){
    est_model <- cICS_mv(data, grid_use, niter, nburn, m0, k0, S0, n0, mass, napprox, nupd,
                         out_param, out_dens, sigma_PY, print_message, light_dens)
  } else if(method == "SLI"){
    est_model <- cSLI_mv(data, grid_use, niter, nburn, m0, k0, S0, n0,
                         mass, nupd, out_param, out_dens, sigma_PY, print_message, light_dens)
  } else if(method == "MAR"){
    est_model <- MAR_mv(data, grid_use, niter, nburn, m0, k0, S0, n0,
                        mass, nupd, out_param, out_dens, sigma_PY, print_message, light_dens)
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
                        univariate = FALSE)
    }else{
      output <- modCond(clust = est_model$clust,
                        niter = niter,
                        nburn = nburn,
                        nnew = as.vector(est_model$newval),
                        tot_time = est_model$time,
                        univariate = FALSE)
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
                        univariate = FALSE)
    }else{
      output <- modCond(clust = est_model$clust,
                        mean = est_model$mu,
                        sigma2 = est_model$s2,
                        probs = est_model$probs,
                        niter = niter,
                        nburn = nburn,
                        nnew = as.vector(est_model$newval),
                        tot_time = est_model$time,
                        univariate = FALSE)
    }
  }

  return(output)
}
