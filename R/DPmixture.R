# options(add.error.underscore=FALSE)
#' Estimating the posterion mean of a Dirichlet Process Mixture Model via marginal Gibbs sampler.
#' The base measure is specified with independent component.
#'
#' @param data A multivariate dataset.
#' @param grid An expanded grid of points to evaluating the density.
#' @param MCMC_param A list of MCMC parameters
#' @inheritParams nsim Number of iterations (includes the burnin iteration)
#' @inheritParams nburn Number of burn-in iterations
#' @inheritParams napprox Number of simulated values to approximate the predictive distribution when updating the latent components
#'
#' @return A list of
#'

DPmixMulti <- function(data = NULL,
                       grid = NULL,
                       MCMC_par = list(),
                       starting_val = list(),
                       params = list(),
                       fix = FALSE,
                       seed = 42) {

  data = as.matrix(data)
  grid = as.matrix(grid)

  #-----------------#
  # MCMC parameters #
  #-----------------#

  nsim = MCMC_par$nsim
  nburn = MCMC_par$nburn
  napprox = MCMC_par$napprox
  nupd = MCMC_par$nupd

  #-----------------#
  # starting values #
  #-----------------#

  conf_start = starting_val$conf_start
  mu_start = starting_val$mu_start
  Lambda_start = starting_val$Lambda_start
  theta_start = starting_val$theta_start
  m0 = starting_val$m0
  B0 = starting_val$B0
  sigma = starting_val$sigma
  theta_fix = starting_val$theta_fix

  #-------#
  # param #
  #-------#

  nu0 = params$nu0
  b1 = params$b1
  B1 = params$B1
  m1 = params$m1
  M1 = params$M1
  s1 = params$s1
  S1 = params$S1
  t1 = params$t1
  t2 = params$t2

  #--------------------#
  # seed and dimension #
  #--------------------#

  set.seed(seed = seed)
  d = ncol(data)
  grid_l = nrow(grid)

  #---------------------------------#
  # check on parameters: STOP check #
  #---------------------------------#

  if(is.null(nsim) || is.null(nburn)) stop("One or more MCMC parameters are missing")
  if(is.null(data)) stop("Give me a data set!")
  if(is.null(grid)) stop("Give me a grid!")
  if(isTRUE(fix) && is.null(theta_fix)) stop("Missing value for theta in params")

  #---------------------------------------------#
  # check on parameters: warning and initialize #
  #---------------------------------------------#

  if(is.null(napprox)) napprox = 100
  if(is.null(nupd)) nupd = round(nsim * .1)
  if(is.null(conf_start)) conf_start = rep(0, nrow(data))
  if(is.null(mu_start)) mu_start = colMeans(data)
  if(is.null(Lambda_start)) Lambda_start = var(data)
  if(is.null(theta_start) && isTRUE(!fix)) theta_start = 1
  if(is.null(theta_start) && isTRUE(fix)) theta_start = theta_fix
  if(is.null(m0)) m0 = colMeans(data)
  if(is.null(B0)) B0 = var(data)
  if(is.null(sigma)) sigma = var(data)
  if(is.null(t1)) t1 = 1
  if(is.null(t2)) t2 = 1
  if(is.null(nu0)) nu0 = d + 2
  if(is.null(b1)) b1 = d + 2
  if(is.null(s1)) s1 = d + 2
  if(is.null(B1)) B1 = diag(1,d)
  if(is.null(M1)) M1 = diag(1,d)
  if(is.null(S1)) S1 = diag(1,d)
  if(is.null(m1)) m1 = rep(0,d)

  if(is.null(m0) || is.null(B0) || is.null(sigma) || is.null(t1) || is.null(t2) || is.null(nu0) ||
     is.null(b1) || is.null(s1) || is.null(B1) || is.null(M1) || is.null(S1) || is.null(m1)){
    warning("One or more parameter missed, initialized to default.")
  }

  #----------------------#
  # Estimating the model #
  #----------------------#

  mod <- marginal_DP_multi_indep(nsim = nsim,
                                 nburn = nburn,
                                 napprox = napprox,
                                 d = d,
                                 data = data,
                                 grid = grid,
                                 conf_start = conf_start,
                                 mu_start = mu_start,
                                 Lambda_start = Lambda_start,
                                 theta = theta_start,
                                 m0 = m0,
                                 B0 = B0,
                                 nu0 = nu0,
                                 sigma = sigma,
                                 b1 = b1,
                                 B1 = B1,
                                 m1 = m1,
                                 M1 = M1,
                                 s1 = s1,
                                 S1 = S1,
                                 t1 = t1,
                                 t2 = t2,
                                 nupd = nupd,
                                 FIX = fix)
  return(mod)
}
