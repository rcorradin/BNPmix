#' @name condMCMC
#' @export condMCMC
#' 
#' @title Estimate an univariate Dirichlet process mixture model with Gaussian kernel using the 
#' conditional Polya urn scheme or Slice sampler
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
#' @param method Different methods to estimate the model. Possible vaule "CPU" or "SLI", conditional Polya urn scheme or Slice sampler, default "CPU".
#' @param napprox Number of approximating value for the conditioning distribution via "CPU" method, default 100.
#' @param nupd How frequently show the curren state of the estimation (number of iterations) - default 1000.
#' @param out_param If TRUE, save the parameters for each iteration, default FALSE.
#' @param out_dens If TRUE, save the parameters for each iteration, default TRUE
#'
#' @return A modCond class object contain the estimated density for each iterations, the allocations for each iterations. If out_param is TRUE, also the parameters.
#' 
#' @examples 
#' data_toy <- c(rnorm(100, -3, 1), rnorm(100, 3, 1))
#' grid <- seq(-7, 7, length.out = 50)
#' est_model <- condMCMC(data = data_toy, grid = grid, niter = 5000, 
#'                       nburn = 1000, napprox = 100)
#' plot(est_model)
#' 

condMCMC <- function(data, grid = NULL, niter, nburn, m0 = NULL, k0 = NULL, 
                     a0 = NULL, b0 = NULL, mass = 1, method = "CPU", 
                     napprox = 100, nupd = 1000, out_param = F, out_dens= TRUE,
                     process = "DP", sigma_PY = 0, print_message = TRUE){
  
  if(is.null(grid)) grid <- 0
  grid_use <- as.vector(grid)
  if(length(dim(data)) > 1){
    stop("The dataset must be a vector")
  }
  if(is.null(m0)) m0 <- 0
  if(is.null(k0)) k0 <- 1
  if(is.null(a0)) a0 <- 2
  if(is.null(b0)) b0 <- 1
  
  if(method == "CPU"){
    if(process == "DP"){
      est_model <- cPUS(data, grid_use, niter, nburn, m0, k0, a0, b0, mass, 
                        napprox, nupd, out_param, out_dens, 0, sigma_PY, print_message)
    } else if(process == "PY"){
      est_model <- cPUS(data, grid_use, niter, nburn, m0, k0, a0, b0, mass, 
                        napprox, nupd, out_param, out_dens, 1, sigma_PY, print_message)
    } 
  } else if(method == "SLI"){
    if(process == "DP"){
      est_model <- cSLI(data, grid_use, niter, nburn, m0, k0, a0, b0, 
                        mass, nupd, out_param, out_dens, 0, sigma_PY, print_message)  
    } else if(process == "PY"){
      est_model <- cSLI(data, grid_use, niter, nburn, m0, k0, a0, b0, 
                        mass, nupd, out_param, out_dens, 1, sigma_PY, print_message)  
    }
  } else if(method == "MPU"){
    if(process == "DP"){
      est_model <- MPU(data, grid_use, niter, nburn, m0, k0, a0, b0, 
                       mass, nupd, out_param, out_dens, 0, sigma_PY, print_message)  
    } else if(process == "PY"){
      est_model <- MPU(data, grid_use, niter, nburn, m0, k0, a0, b0, 
                        mass, nupd, out_param, out_dens, 1, sigma_PY, print_message)  
    }
  }
  
  if(!isTRUE(out_param)){
    if(isTRUE(out_dens)){
      output <- new(Class = "modCond", 
                    density = est_model$dens, 
                    grideval = grid_use, 
                    clust = est_model$clust, 
                    niter = niter, 
                    nburn = nburn, 
                    nnew = as.vector(est_model$newval), 
                    nclust = as.vector(est_model$nclust), 
                    tot_time = est_model$time)
    }else{
      output <- new(Class = "modCond", 
                    clust = est_model$clust, 
                    niter = niter, 
                    nburn = nburn, 
                    nnew = as.vector(est_model$newval), 
                    nclust = as.vector(est_model$nclust), 
                    tot_time = est_model$time)
    }
  } else {
    if(isTRUE(out_dens)){
      output <- new(Class = "modCond", 
                    density = est_model$dens, 
                    grideval = grid_use, 
                    clust = est_model$clust, 
                    mean = est_model$mu, 
                    sigma2 = est_model$s2, 
                    probs = est_model$probs, 
                    niter = niter, 
                    nburn = nburn, 
                    nnew = as.vector(est_model$newval), 
                    nclust = as.vector(est_model$nclust),  
                    tot_time = est_model$time)
    }else{
      output <- new(Class = "modCond", 
                    clust = est_model$clust, 
                    mean = est_model$mu, 
                    sigma2 = est_model$s2, 
                    probs = est_model$probs, 
                    niter = niter, 
                    nburn = nburn, 
                    nnew = as.vector(est_model$newval), 
                    nclust = as.vector(est_model$nclust),  
                    tot_time = est_model$time)
    }
  }
  
  return(output)
}

# METHODS for class modCond ------------------------------------------------------------------------------

setMethod(f = "plot",
          signature(x = "modCond"),
          definition = function(x, prob = c(0.025, 0.975), col = "#0037c4"){
            stopifnot(class(x) == "modCond")
            
            plot_df <- as.data.frame(cbind(x@grideval, colMeans(x@density), 
                                           apply(x@density, 2, function(y) quantile(y, probs = prob[1])),
                                           apply(x@density, 2, function(y) quantile(y, probs = prob[2]))))
            names(plot_df) = c("V1", "V2", "V3", "V4")
            ggplot2::ggplot(plot_df, aes(x = V1, y = V2)) + 
              theme_bw() +
              theme(axis.ticks = element_blank(),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank()) + 
              geom_ribbon(aes(ymin = V3, ymax = V4), alpha = 0.3, fill = col) +
              geom_line(aes(x = V1, y = V2), size= 1, color = col)
            
              # geom_line(aes(x = V1, y = V2), size=.5, col = col) + 
              # geom_line(aes(x = V1, y = V3), size=.5, linetype="dashed", col = col) +
              # geom_line(aes(x = V1, y = V4), size=.5, linetype="dashed", col = col) 
          })


setMethod(f = "acfplot",
          signature(x = "modCond"),
          definition = function(x, maxlag = 30, alpha = 0.95){
            stopifnot(class(x) == "modCond")
            
            acf_temp <- acf(x = apply(x@clust, 1, function(y) length(unique(y))), lag.max = maxlag, plot = FALSE)
            conf.lims <- c(-1,1)*qnorm((1 + alpha)/2)/sqrt(acf_temp$n.used)
            plot_db <- data.frame(acf_temp$lag, acf_temp$acf)
            ggplot2::ggplot(plot_db, aes(x=acf_temp.lag, y = acf_temp.acf)) + scale_x_continuous(breaks=seq(0,maxlag + 1,round(maxlag / 4))) +
              theme_bw() +
              geom_hline(yintercept=conf.lims, lty=2, col='blue') +
              labs(y="Autocorrelations", x="Lag", title= "") +
              geom_segment(aes(xend=acf_temp.lag, yend=0)) + geom_point()
            
          })