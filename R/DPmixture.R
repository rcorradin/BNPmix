# options(add.error.underscore=FALSE)
DPmixture <- function(nsim, nburn, napprox, nparam, d, 
                  grid_l, data, grid, mu_start, Lambda_start, theta_start, 
                  m0, B0, nu0, sigma, b1, B1, m1, M1, s1, S1, 
                  t1, t2, nupd, plim, fix = TRUE, seed = 42){
  
  set.seed(seed = seed)
  
  if(fix != TRUE){
    mod <- main_fun(nsim = nsim,
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