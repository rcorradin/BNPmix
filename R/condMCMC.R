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
#' @param out_dens If TRUE, return also the estimated density, default TRUE
#' @param process Dirichlet process ("DP") or Pitman-Yor process ("PY"), default "DP"
#' @param sigma_PY Discount parameter of the Pitman-Yor process, default 0
#' @param print_message If TRUE print the status of the estimation, default TRUE
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
                     a0 = NULL, b0 = NULL, mass = 1, method = "ICS",
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

  if(method == "ICS"){
    if(process == "DP"){
      est_model <- cICS(data, grid_use, niter, nburn, m0, k0, a0, b0, mass,
                        napprox, nupd, out_param, out_dens, 0, sigma_PY, print_message)
    } else if(process == "PY"){
      est_model <- cICS(data, grid_use, niter, nburn, m0, k0, a0, b0, mass,
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
  } else if(method == "MAR"){
    if(process == "DP"){
      est_model <- MAR(data, grid_use, niter, nburn, m0, k0, a0, b0,
                       mass, nupd, out_param, out_dens, 0, sigma_PY, print_message)
    } else if(process == "PY"){
      est_model <- MAR(data, grid_use, niter, nburn, m0, k0, a0, b0,
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
                    tot_time = est_model$time)
    }else{
      output <- new(Class = "modCond",
                    clust = est_model$clust,
                    niter = niter,
                    nburn = nburn,
                    nnew = as.vector(est_model$newval),
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
            ggplot2::ggplot(plot_df, mapping = ggplot2::aes(x = V1, y = V2)) +
              ggplot2::theme_bw() +
              ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                              axis.title.x = ggplot2::element_blank(),
                              axis.title.y = ggplot2::element_blank()) +
              ggplot2::geom_ribbon(mapping = ggplot2::aes(ymin = V3, ymax = V4), alpha = 0.3, fill = col) +
              ggplot2::geom_line(mapping = ggplot2::aes(x = V1, y = V2), size= 1, color = col)
          })
