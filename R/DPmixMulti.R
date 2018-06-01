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
#' @inheritSection estimated density
#' @inheritSection partition matrix
#' @inheritSection mass values
#'

DPmixMulti <- function(nsim = NULL,
                       nburn = NULL,
                       napprox = 100,
                       data = NULL,
                       grid = NULL,
                       conf_start = NULL,
                       mu_start = NULL,
                       Lambda_start = NULL,
                       theta = NULL,
                       m0 = NULL,
                       B0 = NULL,
                       nu0 = NULL,
                       sigma = NULL,
                       b1 = NULL,
                       B1 = NULL,
                       m1 = NULL,
                       M1 = NULL,
                       s1 = NULL,
                       S1 = NULL,
                       t1 = NULL,
                       t2 = NULL,
                       k1 = NULL,
                       nupd = NULL,
                       fix = FALSE,
                       dep = FALSE,
                       seed = 42) {

  if(is.null(data)) stop("Give me a data set!")
  if(is.null(grid)) stop("Give me a grid!")

  data = as.matrix(data)
  grid = as.matrix(grid)

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
  if(isTRUE(fix) && is.null(theta)) stop("Missing value for theta in params")

  #---------------------------------------------#
  # check on parameters: warning and initialize #
  #---------------------------------------------#

  if(is.null(m0) || is.null(B0) || is.null(sigma) || (is.null(t1) & isTRUE(!fix)) ||
     (is.null(t2) & isTRUE(!fix)) || is.null(nu0) ||
     is.null(b1) || (is.null(s1) & isTRUE(!dep)) || is.null(B1) || is.null(M1) ||
     (is.null(S1) & isTRUE(!dep)) || is.null(m1) || (is.null(k1) & isTRUE(dep))){
    warning("One or more parameter missed, initialized to default.")
  }

  if(is.null(nupd)) nupd = round(nsim * .1)
  if(is.null(conf_start)) conf_start = rep(0, nrow(data))
  if(is.null(mu_start)) mu_start = colMeans(data)
  if(is.null(Lambda_start)) Lambda_start = var(data)
  if(is.null(theta) && isTRUE(!fix)) theta = 1
  if(is.null(theta) && isTRUE(fix)) theta = theta
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
  if(is.null(k1)) k1 = 1

  #----------------------#
  # Estimating the model #
  #----------------------#

  if(dep == FALSE){
    mod <- marginal_DP_multi_indep_indep(nsim = nsim,
                                         nburn = nburn,
                                         napprox = napprox,
                                         d = d,
                                         data = data,
                                         grid = grid,
                                         conf_start = conf_start,
                                         mu_start = mu_start,
                                         Lambda_start = Lambda_start,
                                         theta = theta,
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

  }else{
    mod <- marginal_DP_multi_indep_dep(nsim = nsim,
                                      nburn = nburn,
                                      napprox = napprox,
                                      d = d,
                                      data = data,
                                      grid = grid,
                                      conf_start = conf_start,
                                      mu_start = mu_start,
                                      Lambda_start = Lambda_start,
                                      theta = theta,
                                      m0 = m0,
                                      B0 = B0,
                                      nu0 = nu0,
                                      sigma = sigma,
                                      b1 = b1,
                                      B1 = B1,
                                      m1 = m1,
                                      k1 = k1,
                                      t1 = t1,
                                      t2 = t2,
                                      nupd = nupd,
                                      FIX = fix)
  }
  return(mod)
}
