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
#' @param out_dens If TRUE, save the parameters for each iteration, default TRUE
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
                       S0 = NULL, n0 = NULL, mass = 1, method = "CPU",
                       napprox = 100, nupd = 1000, out_param = F, out_dens = TRUE,
                       process = "DP", sigma_PY = 0, print_message = TRUE, light_dens = TRUE){

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

  if(method == "CPU"){
    if(process == "DP"){
      est_model <- cPUS_mv(data, grid_use, niter, nburn, m0, k0, S0, n0,
                           mass, napprox, nupd, out_param, out_dens, 0, sigma_PY, print_message, light_dens)
    } else if(process == "PY"){
      est_model <- cPUS_mv(data, grid_use, niter, nburn, m0, k0, S0, n0,
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
  } else if(method == "MPU"){
    if(process == "DP"){
      est_model <- MPU_mv(data, grid_use, niter, nburn, m0, k0, S0, n0,
                           mass, nupd, out_param, out_dens, 0, sigma_PY, print_message, light_dens)
    } else if(process == "PY"){
      est_model <- MPU_mv(data, grid_use, niter, nburn, m0, k0, S0, n0,
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
                    nclust = as.vector(est_model$nclust),
                    tot_time = est_model$time)
    }else{
      output <- new(Class = "modCondMv",
                    clust = est_model$clust,
                    niter = niter,
                    nburn = nburn,
                    nnew = as.vector(est_model$newval),
                    nclust = as.vector(est_model$nclust),
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
                    nclust = as.vector(est_model$nclust),
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
                    nclust = as.vector(est_model$nclust),
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
                ggplot2::ggplot(plot_df_use, aes(x = Group.1, y = V1)) +
                  theme_bw() +
                  theme(axis.ticks = element_blank(),
                        axis.title.x = element_blank(),
                        axis.title.y = element_blank()) +
                  geom_ribbon(aes(ymin = V2, ymax = V3), alpha = 0.3, fill = col) +
                  geom_line(aes(x = Group.1, y = V1), size= 1, color = col)

                # geom_line(aes(x = Group.1, y = V1), size=.5, col = col) +
                # geom_line(aes(x = Group.1, y = V2), size=.5, linetype="dashed", col = col) +
                # geom_line(aes(x = Group.1, y = V3), size=.5, linetype="dashed", col = col)
              }else{
                plot_df_use <- aggregate(plot_df, by = list(plot_df[[dimension[1]]],plot_df[[dimension[2]]]), FUN = sum)
                ggplot2::ggplot(data = plot_df_use, mapping = aes(x = Group.1, y = Group.2, z = V1)) +
                  stat_contour(data = plot_df_use, mapping = aes(x = Group.1, y = Group.2, z = V1), bins = 10, col = col) +
                  theme_bw() +
                  theme(axis.ticks = element_blank(),
                        axis.title.x = element_blank(),
                        axis.title.y = element_blank())
              }
            } else {
              plot_df <- as.data.frame(cbind(x@grideval, x@density))
              names(plot_df) = c(paste("GR", 1:ncol(x@grideval), sep = ''), "V1")

              if(dimension[1] == dimension[2]){
                plot_df_use <- aggregate(plot_df, by = list(plot_df[[dimension[1]]]), FUN = sum)
                ggplot2::ggplot(data = plot_df_use, aes(x = Group.1, y = V1)) +
                  theme_bw() +
                  theme(axis.ticks = element_blank(),
                        axis.title.x = element_blank(),
                        axis.title.y = element_blank()) +
                  geom_line(aes(x = Group.1, y = V1), size= 1, color = col)

                # geom_line(aes(x = Group.1, y = V1), size=.5, col = col) +
                # geom_line(aes(x = Group.1, y = V2), size=.5, linetype="dashed", col = col) +
                # geom_line(aes(x = Group.1, y = V3), size=.5, linetype="dashed", col = col)
              }else{
                plot_df_use <- aggregate(plot_df, by = list(plot_df[[dimension[1]]],plot_df[[dimension[2]]]), FUN = sum)
                ggplot2::ggplot(data = plot_df_use, mapping = aes(x = Group.1, y = Group.2, z = V1)) +
                  stat_contour(data = plot_df_use, mapping = aes(x = Group.1, y = Group.2, z = V1), bins = 10, col = col) +
                  theme_bw() +
                  theme(axis.ticks = element_blank(),
                        axis.title.x = element_blank(),
                        axis.title.y = element_blank())
              }
            }
          })

setMethod(f = "acfplot",
          signature(x = "modCondMv"),
          definition = function(x, maxlag = 30, alpha = 0.95){
            stopifnot(class(x) == "modCondMv")

            acf_temp <- acf(x = apply(x@clust, 1, function(y) length(unique(y))), lag.max = maxlag, plot = FALSE)
            conf.lims <- c(-1,1)*qnorm((1 + alpha)/2)/sqrt(acf_temp$n.used)
            plot_db <- data.frame(acf_temp$lag, acf_temp$acf)
            ggplot2::ggplot(plot_db, aes(x=acf_temp.lag, y = acf_temp.acf)) + scale_x_continuous(breaks=seq(0,maxlag + 1,round(maxlag / 4))) +
                            theme_bw() +
                            geom_hline(yintercept=conf.lims, lty=2, col='blue') +
                            labs(y="Autocorrelations", x="Lag", title= "") +
                            geom_segment(aes(xend=acf_temp.lag, yend=0)) + geom_point()

          })
