# options(add.error.underscore=FALSE)
DPmixture <- function(nsim = NULL,
                      nburn = NULL,
                      napprox = NULL,
                      nparam = NULL,
                      data = NULL,
                      grid = NULL,
                      mu_start = NULL,
                      Lambda_start = NULL,
                      theta_start = NULL,
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
                      nupd = NULL,
                      plim = NULL,
                      fix = TRUE,
                      seed = 42) {

  set.seed(seed = seed)
  d = ncol(data)
  grid_l = nrow(grid)

  # check if any parameter is missed

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

  if(is.null(data)){
    stop("Give me a data set!")
  }

  if(is.null(grid)){
    stop("Give me a grid!")
  }

  if(is.null(mu_start)){
    mu_start = colMeans(data)
  }

  if(is.null(Lambda_start)){
    Lambda_start = var(data)
  }



  if(fix != TRUE){
    mod <- main_fun(nsim = nsim,
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
                    plim = plim)
  }


  if(fix == TRUE){
    mod <- main_fun_fix(nsim = nsim,
                    nburn = nburn,
                    napprox = napprox,
                    nparam = nparam,
                    d = d,
                    grid_l = grid_l,
                    data = data,
                    grid = grid,
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
                    nupd = nupd,
                    plim = plim)
  }

  return(mod)
}
