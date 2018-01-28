# options(add.error.underscore=FALSE)
DPmixMulti <- function(data = NULL,
                       grid = NULL,
                       MCMC_par = list(nsim= NULL, nburn = NULL, napprox = 100, nparam, nupd, plim),
                       starting_val = list(mu_start = NULL, Lambda_start = NULL, theta_start = NULL,
                                           m0 = NULL, B0 = NULL, sigma = NULL),
                       params = list(nu0 = NULL, b1 = NULL, B1 = NULL, m1 = NULL, M1 = NULL,
                                     s1 = NULL, S1 = NULL, t1 = NULL, t2 = NULL, theta_fix = NULL),
                       fix = FALSE,
                       seed = 42) {

  #-----------------#
  # MCMC parameters #
  #-----------------#

  nsim= MCMC_par$nsim
  nburn = MCMC_par$nburn
  napprox = MCMC_par$napprox
  nparam = MCMC_par$nparam
  nupd = MCMC_par$nupd
  plim = MCMC_par$plim

  #-----------------#
  # starting values #
  #-----------------#

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
  t2 = params$r2

  #--------------------#
  # seed and dimension #
  #--------------------#

  set.seed(seed = seed)
  d = ncol(data)
  grid_l = nrow(grid)

  #--------------------------------#
  # check on parameters: MCMCparam #
  #--------------------------------#

  if(is.null(nsim) || is.null(nburn)){
    stop("One or more MCMC parameters are missing")
  }

  if(is.null(napprox) || is.null(nparam) || is.null(plim)) {
    warning("napprox or nparam or plim missed, using the default one")
    if(is.null(napprox)){
      napprox = 100
    }
    if(is.null(nparam)){
      nparam = round(nrow(data) * 1.5)
    }
    if(is.null(plim)){
      plim = round(0.1 * nrow(data))
    }
  }

  if(is.null(nupd)){
    nupd = round(nsim * .1)
  }

  #---------------------------#
  # check on parameters: data #
  #---------------------------#

  if(is.null(data)){
    stop("Give me a data set!")
  }

  if(is.null(grid)){
    stop("Give me a grid!")
  }

  #-----------------------------------#
  # check on parameters: starting_val #
  #-----------------------------------#

  if(is.null(mu_start)){
    mu_start = colMeans(data)
  }

  if(is.null(Lambda_start)){
    Lambda_start = var(data)
  }

  if(is.null(theta_start)){
    if(isTRUE(fix)){
      theta_start = theta_fix
    }
    if(isTRUE(!fix)){
      theta_start = 1
      warning("Initialized theta_start = 1")
    }
  }

  if(is.null(m0)){
    m0 = colMeans(data)
  }

  if(is.null(B0)){
    B0 = var(data)
  }

  if(is.null(sigma)){
    sigma = var(data)
  }

  #-----------------------------#
  # check on parameters: params #
  #-----------------------------#

  if(is.null(t1) || is.null(t2)){
    t1 = 1
    t2 = 1
    if(isTRUE(!fix)){
      warning("Parameters t1 and t2 missed: initialized to default (equal to 1)")
    }
  }

  if(isTRUE(fix) && is.null(theta_fix)){
    stop("Missing value for theta in params")
  }

  if(is.null(nu0)){
    nu0 = d + 2
    warning("Missing nu0 in params, initialized to minimum (ncol(data)+2)")
  }

  if(is.null(b1)){
    b1 = d + 2
    warning("Missing b1 in params, initialized to minimum (ncol(data)+2)")
  }

  if(is.null(s1)){
    s1 = d + 2
    warning("Missing s1 in params, initialized to minimum (ncol(data)+2)")
  }

  if(is.null(B1)){
    B1 = diag(1,d)
    warning("Missing B1 in params, initialized to default (diag(1))")
  }

  if(is.null(M1)){
    M1 = diag(1,d)
    warning("Missing M1 in params, initialized to default (diag(1))")
  }

  if(is.null(S1)){
    S1 = diag(1,d)
    warning("Missing S1 in params, initialized to default (diag(1))")
  }

  if(is.null(m1)){
    m1 = rep(0,d)
    warning("Missing m1 in params, initialized to default (diag(1))")
  }

  #----------------------#
  # Estimating the model #
  #----------------------#

  mod <- marginal_DP_multi(nsim = nsim,
                           nburn = nburn,
                           napprox = napprox,
                           nparam = nparam,
                           d = d,
                           grid_l = grid_l,
                           data = data,
                           grid = grid,
                           conf_start,
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
                           plim = plim,
                           FIX = fix)
  return(mod)
}
