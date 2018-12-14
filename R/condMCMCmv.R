#' @name condMCMCmv
#' @export condMCMCmv
#'
#' @title Estimate a multivariate Dirichlet process mixture model with Gaussian kernel using the
#' conditional Polya urn scheme or Slice sampler
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
#' @param method Different methods to estimate the model. Possible vaule "CPU" or "SLI", conditional Polya urn scheme or Slice sampler, default "CPU".
#' @param napprox Number of approximating value for the conditioning distribution via "CPU" method, default 100.
#' @param nupd How frequently show the curren state of the estimation (number of iterations) - default 1000.
#' @param out_param If TRUE, save the parameters for each iteration, default FALSE.
#' @param out_dens If TRUE, return also the estimated density, default TRUE
#' @param process Dirichlet process ("DP") or Pitman-Yor process ("PY"), default "DP"
#' @param sigma_PY Discount parameter of the Pitman-Yor process, default 0
#' @param print_message If TRUE print the status of the estimation, default TRUE
#' @param light_dens If TRUE return only the mean of the estimated densities, default TRUE
#'
#' @return A modCond class object contain the estimated density for each iterations, the allocations for each iterations. If out_param is TRUE, also the parameters.
#'
#' @examples
#' data_toy <- cbind(c(rnorm(100, -3, 1), rnorm(100, 3, 1)),
#'                   c(rnorm(100, -3, 1), rnorm(100, 3, 1)))
#' grid <- expand.grid(seq(-7, 7, length.out = 50),
#'                     seq(-7, 7, length.out = 50))
#' est_model <- condMCMCmv(data = data_toy, grid = grid, niter = 5000,
#'                        nburn = 1000, napprox = 100)
#' plot(est_model)
#'

condMCMCmv <- function(data, grid = NULL, niter, nburn,  m0 = NULL, k0 = NULL,
                       S0 = NULL, n0 = NULL, mass = 1, method = "ICS",
                       napprox = 100, nupd = 1000, out_param = F, out_dens = TRUE,
                       process = "DP", sigma_PY = 0, print_message = TRUE,
                       light_dens = TRUE){

  if(is.null(grid)) grid <- matrix(0, ncol = ncol(data))
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
    if(process == "DP"){
      est_model <- cICS_mv(data, grid_use, niter, nburn, m0, k0, S0, n0,
                           mass, napprox, nupd, out_param, out_dens, 0, sigma_PY, print_message, light_dens)
    } else if(process == "PY"){
      est_model <- cICS_mv(data, grid_use, niter, nburn, m0, k0, S0, n0,
                           mass, napprox, nupd, out_param, out_dens, 1, sigma_PY, print_message, light_dens)
    }
  } else if(method == "SLI"){
    if(process == "DP"){
      est_model <- cSLI_mv(data, grid_use, niter, nburn, m0, k0, S0, n0,
                           mass, nupd, out_param, out_dens, 0, sigma_PY, print_message, light_dens)
    } else if(process == "PY"){
      est_model <- cSLI_mv(data, grid_use, niter, nburn, m0, k0, S0, n0,
                           mass, nupd, out_param, out_dens, 1, sigma_PY, print_message, light_dens)
    }
  } else if(method == "MAR"){
    if(process == "DP"){
      est_model <- MAR_mv(data, grid_use, niter, nburn, m0, k0, S0, n0,
                           mass, nupd, out_param, out_dens, 0, sigma_PY, print_message, light_dens)
    } else if(process == "PY"){
      est_model <- MAR_mv(data, grid_use, niter, nburn, m0, k0, S0, n0,
                           mass, nupd, out_param, out_dens, 1, sigma_PY, print_message, light_dens)
    }
  }

  if(!isTRUE(out_param)){
    if(isTRUE(out_dens)){
      output <- new(Class = "modCondMv",
                    density = est_model$dens,
                    grideval = grid_use,
                    clust = est_model$clust,
                    niter = niter,
                    nburn = nburn,
                    nnew = as.vector(est_model$newval),
                    tot_time = est_model$time)
    }else{
      output <- new(Class = "modCondMv",
                    clust = est_model$clust,
                    niter = niter,
                    nburn = nburn,
                    nnew = as.vector(est_model$newval),
                    tot_time = est_model$time)
    }
  } else {
    if(isTRUE(out_dens)){
      output <- new(Class = "modCondMv",
                    density = est_model$dens,
                    grideval = grid_use,
                    clust = est_model$clust,
                    mean = as.list(est_model$mu),
                    sigma2 = est_model$s2,
                    probs = est_model$probs,
                    niter = niter,
                    nburn = nburn,
                    nnew = as.vector(est_model$newval),
                    tot_time = est_model$time)
    }else{
      output <- new(Class = "modCondMv",
                    clust = est_model$clust,
                    mean = as.list(est_model$mu),
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

# METHODS for class modCondMv ----------------------------------------------------------------------------

setMethod(f = "plot",
          signature(x = "modCondMv"),
          definition = function(x, prob = c(0.025, 0.975), dimension = c(1,2), col = "#0037c4"){
            stopifnot(class(x) == "modCondMv")

            if(dim(x@density)[2] > 1){
              plot_df <- as.data.frame(cbind(x@grideval, colMeans(x@density),
                                             apply(x@density, 2, function(y) quantile(y, probs = prob[1])),
                                             apply(x@density, 2, function(y) quantile(y, probs = prob[2]))))
              names(plot_df) = c(paste("GR", 1:ncol(x@grideval), sep = ''), "V1", "V2", "V3")

              if(dimension[1] == dimension[2]){
                plot_df_use <- aggregate(plot_df, by = list(plot_df[[dimension[1]]]), FUN = sum)
                ggplot2::ggplot(plot_df_use, mapping = ggplot2::aes(x = Group.1, y = V1)) +
                  ggplot2::theme_bw() +
                  ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                                  axis.title.x = ggplot2::element_blank(),
                                  axis.title.y = ggplot2::element_blank()) +
                  ggplot2::geom_ribbon(mapping = ggplot2::aes(ymin = V2, ymax = V3), alpha = 0.3, fill = col) +
                  ggplot2::geom_line(mapping = ggplot2::aes(x = Group.1, y = V1), size= 1, color = col)

              }else{
                plot_df_use <- aggregate(plot_df, by = list(plot_df[[dimension[1]]],plot_df[[dimension[2]]]), FUN = sum)
                ggplot2::ggplot(data = plot_df_use, mapping = ggplot2::aes(x = Group.1, y = Group.2, z = V1)) +
                  ggplot2::stat_contour(data = plot_df_use, mapping = ggplot2::aes(x = Group.1, y = Group.2, z = V1), bins = 10, col = col) +
                  ggplot2::theme_bw() +
                  ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                                  axis.title.x = ggplot2::element_blank(),
                                  axis.title.y = ggplot2::element_blank())
              }
            } else {
              plot_df <- as.data.frame(cbind(x@grideval, x@density))
              names(plot_df) = c(paste("GR", 1:ncol(x@grideval), sep = ''), "V1")

              if(dimension[1] == dimension[2]){
                plot_df_use <- aggregate(plot_df, by = list(plot_df[[dimension[1]]]), FUN = sum)
                ggplot2::ggplot(data = plot_df_use, mapping = ggplot2::aes(x = Group.1, y = V1)) +
                  ggplot2::theme_bw() +
                  ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                        axis.title.x = ggplot2::element_blank(),
                        axis.title.y = ggplot2::element_blank()) +
                  ggplot2::geom_line(mapping = ggplot2::aes(x = Group.1, y = V1), size= 1, color = col)

                # geom_line(aes(x = Group.1, y = V1), size=.5, col = col) +
                # geom_line(aes(x = Group.1, y = V2), size=.5, linetype="dashed", col = col) +
                # geom_line(aes(x = Group.1, y = V3), size=.5, linetype="dashed", col = col)
              }else{
                plot_df_use <- aggregate(plot_df, by = list(plot_df[[dimension[1]]],plot_df[[dimension[2]]]), FUN = sum)
                ggplot2::ggplot(data = plot_df_use, mapping = ggplot2::aes(x = Group.1, y = Group.2, z = V1)) +
                  ggplot2::stat_contour(data = plot_df_use, mapping = ggplot2::aes(x = Group.1, y = Group.2, z = V1), bins = 10, col = col) +
                  ggplot2::theme_bw() +
                  ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                        axis.title.x = ggplot2::element_blank(),
                        axis.title.y = ggplot2::element_blank())
              }
            }
          })
