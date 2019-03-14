#' @name condDDP
#' @export condDDP
#'
#' @title Conditional dependent Dirichlet process
#' @title Estimate an univariate dependent Dirichlet process mixture model with Gaussian kernel using
#' importance conditional sampler scheme.
#'
#' @param data A dataset (vector).
#' @param group group for the observed data, same dimension as data.
#' @param grid A grid to evaluate the estimated density (vector).
#' @param niter Number of iterations to estimate the model.
#' @param nburn Number of burn-in iterations.
#' @param m0 Mean of the distribution of location component of the base measure.
#' @param k0 Tuning parameter of the location component variance of the base measure.
#' @param a0 Shap parameter of the scale component distribution of the base measure.
#' @param b0 Rate parameter of the scale component distribution of the base measure.
#' @param mass Mass parameter of the Dirichlet process.
#' @param wei weigth of the group process.
#' @param napprox Number of approximating value for the conditioning distribution via "CPU" method, default 100.
#' @param n_approx_unif number of values used in the importance sampling step for the multivariate beta distribution.
#' @param nupd How frequently show the curren state of the estimation (number of iterations) - default 1000.
#' @param out_dens If TRUE, save the parameters for each iteration, default TRUE.
#' @param print_message Print the status of the estimation.
#' @param light_dens Return only the posterior mean of the densities.
#'
#' @return A modCond class object contain the estimated density for each iterations,
#' the allocations for each iterations. If out_param is TRUE, also the parameters.
#'
#' @examples
#' set.seed(42)
#' data_toy <- c(rnorm(50, -4, 1), rnorm(100, 0, 1), rnorm(50, 4, 1))
#' group_toy <- c(rep(1,100), rep(2,100))
#' grid <- seq(-7, 7, length.out = 50)
#' est_model <- condDDP(data = data_toy, group = group_toy, grid = grid, niter = 1000,
#'                      nburn = 100, napprox = 100, nupd = 100)
#' plot(est_model)
#'

condDDP <- function(data, group, grid = NULL, niter, nburn, m0 = NULL, k0 = NULL,
                     a0 = NULL, b0 = NULL, mass = 1, wei = 0.5, napprox = 10,
                     n_approx_unif = 1000, nupd = 1000, out_dens= TRUE,
                    print_message = TRUE, light_dens = FALSE){

  if(is.null(grid)) grid <- 0
  grid_use <- as.vector(grid)
  if(length(dim(data)) > 1){
    stop("The dataset must be a vector")
  }

  group <- as.numeric(as.factor(group))
  ngr   <- length(unique(group))

  if(is.null(m0)) m0 <- 0
  if(is.null(k0)) k0 <- 1
  if(is.null(a0)) a0 <- 2
  if(is.null(b0)) b0 <- 1

  est_model <- cDDP(data,
                    group,
                    ngr,
                    grid,
                    niter,
                    nburn,
                    m0,
                    k0,
                    a0,
                    b0,
                    mass,
                    wei,
                    napprox,
                    n_approx_unif,
                    nupd,
                    out_dens,
                    print_message,
                    light_dens)

  if(isTRUE(out_dens)){
    if(light_dens == TRUE){
      output <- new(Class = "modCondDep",
                    density = est_model$dens[,,1],
                    grideval = grid_use,
                    clust = est_model$clust,
                    group_log = est_model$group_log,
                    niter = niter,
                    nburn = nburn,
                    nclust = as.vector(est_model$nclust),
                    tot_time = est_model$time,
                    group = group,
                    wvals = est_model$wvals)
    } else {
      output <- new(Class = "modCondDep",
                    density = est_model$dens,
                    grideval = grid_use,
                    clust = est_model$clust,
                    group_log = est_model$group_log,
                    niter = niter,
                    nburn = nburn,
                    nclust = as.vector(est_model$nclust),
                    tot_time = est_model$time,
                    group = group,
                    wvals = est_model$wvals)
    }

  }else{
    output <- new(Class = "modCondDep",
                  clust = est_model$clust,
                  group_log = est_model$group_log,
                  niter = niter,
                  nburn = nburn,
                  # nnew = as.vector(est_model$newval),
                  nclust = as.vector(est_model$nclust),
                  tot_time = est_model$time,
                  group = group,
                  wvals = est_model$wvals)
  }

  return(output)
}

# METHODS for class modCond ------------------------------------------------------------------------------

#' @title modCondDep object plot
#' @description \code{plot} method for class \code{modCondDep}
#'
#' @param x object of class \code{modCondDep}.
#' @param ncol number of column in the final plot.

setMethod(f = "plot",
          signature(x = "modCondDep"),
          definition = function(x, ncol = 1){
            with(x, {

              stopifnot(class(x) == "modCondDep")

              ngr <- length(unique(x@group))
              plot_df <- as.data.frame(cbind(rep(x@grideval, ngr), as.vector(apply(x@density, c(1,2), mean)),
                                             as.vector(sapply(1:ngr, function(y) rep(paste("Group ", y), length(x@grideval)))) ))
              plot_df[,1:2] <- as.data.frame(cbind(rep(x@grideval, ngr), as.vector(apply(x@density, c(1,2), mean))))

              ggplot2::ggplot(plot_df, mapping = ggplot2::aes(x = V1, y = V2, color = V3)) +
                ggplot2::theme_bw() +
                ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                               axis.title.x = ggplot2::element_blank(),
                               axis.title.y = ggplot2::element_blank()) +
                ggplot2::geom_line() +
                ggplot2::facet_wrap(~ factor(V3), ncol = 1) +
                ggplot2::guides(fill=FALSE, color=FALSE)
            })
          })

