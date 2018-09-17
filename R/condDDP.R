#' @name condDDP
#' @export condDDP
#' 
#' @title Estimate an univariate Dependent Dirichlet process mixture model with Gaussian kernel using the 
#' conditional Polya urn scheme
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
#' @param out_dens If TRUE, save the parameters for each iteration, default TRUE
#' @param print_message Print the status of the estimation
#'
#' @return A modCond class object contain the estimated density for each iterations, the allocations for each iterations. If out_param is TRUE, also the parameters.
#' 
#' @examples 
#' set.seed(42)
#' data_toy <- c(rnorm(50, -4, 1), rnorm(100, 0, 1), rnorm(50, 4, 1))
#' group_toy <- c(rep(1,100), rep(2,100))
#' grid <- seq(-7, 7, length.out = 50)
#' est_model <- condDDP(data = data_toy, group = group_toy, grid = grid, niter = 10000, 
#'                      nburn = 1000, napprox = 100, a0 = 4)
#' plot(est_model)
#' 

condDDP <- function(data, group, grid = NULL, niter, nburn, m0 = NULL, k0 = NULL, 
                     a0 = NULL, b0 = NULL, mass = 1, wei = 0.5, napprox = 100, 
                     n_approx_unif = 1000, nupd = 1000, out_dens= TRUE, print_message = TRUE){
  
  if(is.null(grid)) grid <- 0
  grid_use <- as.vector(grid)
  if(length(dim(data)) > 1){
    stop("The dataset must be a vector")
  }
  # if(length(data) != length(group)){
  #   stop("The dataset and group must have same length")
  # }
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
                    print_message)
  
  if(isTRUE(out_dens)){
    output <- new(Class = "modCondDep", 
                  density = est_model$dens, 
                  grideval = grid_use, 
                  clust = est_model$clust, 
                  zeta = est_model$zeta, 
                  niter = niter, 
                  nburn = nburn, 
                  nnew = as.vector(est_model$newval), 
                  nclust = as.vector(est_model$nclust), 
                  tot_time = est_model$time,
                  group = group)
  }else{
    output <- new(Class = "modCondDep", 
                  clust = est_model$clust, 
                  zeta = est_model$zeta, 
                  niter = niter, 
                  nburn = nburn, 
                  nnew = as.vector(est_model$newval), 
                  nclust = as.vector(est_model$nclust), 
                  tot_time = est_model$time,
                  group = group)
  }
  
  return(output)
}

# METHODS for class modCond ------------------------------------------------------------------------------

setMethod(f = "plot",
          signature(x = "modCondDep"),
          definition = function(x, prob = c(0.1, 0.9), ncol = 1){
            stopifnot(class(x) == "modCondDep")
            
            ngr <- length(unique(x@group))
            plot_df <- as.data.frame(cbind(rep(x@grideval, ngr), as.vector(apply(x@density, c(1,2), mean)), 
                                           as.vector(apply(x@density, c(1,2), function(y) quantile(y, probs = prob[1]))),
                                           as.vector(apply(x@density, c(1,2), function(y) quantile(y, probs = prob[2]))),
                                           as.vector(sapply(1:ngr, function(y) rep(paste("Group ", y), length(x@grideval)))) ))
            plot_df[,1:4] <- as.data.frame(cbind(rep(x@grideval, ngr), as.vector(apply(x@density, c(1,2), mean)), 
                                           as.vector(apply(x@density, c(1,2), function(y) quantile(y, probs = prob[1]))),
                                           as.vector(apply(x@density, c(1,2), function(y) quantile(y, probs = prob[2])))))
            
            ggplot2::ggplot(plot_df, aes(x = V1, y = V2, color = V5)) + 
              theme_bw() +
              theme(axis.ticks = element_blank(),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank()) + 
              geom_ribbon(aes(ymin = V3, ymax = V4, fill = V5), alpha = 0.3, color = NA) +
              geom_line() +
              facet_wrap(~ factor(V5), ncol = 1) +  guides(fill=FALSE, color=FALSE)
            
            # geom_line(aes(x = V1, y = V2), size=.5, col = col) + 
            # geom_line(aes(x = V1, y = V3), size=.5, linetype="dashed", col = col) +
            # geom_line(aes(x = V1, y = V4), size=.5, linetype="dashed", col = col) 
          })


setMethod(f = "trace_ngr",
          signature(x = "modCondDep"),
          definition = function(x, ncol = 1){
            stopifnot(class(x) == "modCondDep")
            
            ngr <- length(unique(x@group))
            niter <- nrow(x@clust)
            nobs <- ncol(x@clust)
            nclust <- c()
            
            for(i in 0:ngr){
              nclust <- c(nclust, apply(cbind(x@clust, x@zeta), 1, function(y) length(unique(y[1:nobs][y[(nobs+1):(2*nobs)] == i])) ))
            }
            plot_df <- as.data.frame(cbind(rep(1:niter, ngr + 1), nclust, 
                                           as.vector(sapply(0:ngr, function(y) rep(paste("Group ", y), niter)))))
            plot_df[,1:2] <- as.data.frame(cbind(rep(1:niter, ngr + 1), nclust))
            plot_df <- rbind(plot_df[-c(1:niter),], plot_df[c(1:niter),])
            colnames(plot_df) <- c("V1", "V2", "V3")
            ggplot2::ggplot(plot_df, aes(x = V1, y = V2, color = V3)) + 
              theme_bw() +
              theme(axis.ticks = element_blank(),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank()) + 
              geom_line() +
              facet_wrap(~ factor(V3), ncol = 1) +  guides(fill=FALSE, color=FALSE)
            
            # geom_line(aes(x = V1, y = V2), size=.5, col = col) + 
            # geom_line(aes(x = V1, y = V3), size=.5, linetype="dashed", col = col) +
            # geom_line(aes(x = V1, y = V4), size=.5, linetype="dashed", col = col) 
          })

setMethod(f = "trace_obs",
          signature(x = "modCondDep"),
          definition = function(x, ncol = 1){
            stopifnot(class(x) == "modCondDep")
            
            ngr <- length(unique(x@group))
            niter <- nrow(x@clust)
            nobs <- ncol(x@clust)
            nobs_v <- c()
            
            for(i in 0:ngr){
              nobs_v <- c(nobs_v, apply(x@zeta, 1, function(y) sum(y == i)))
            }
            plot_df <- as.data.frame(cbind(rep(1:niter, ngr + 1), nobs_v, 
                                           as.vector(sapply(0:ngr, function(y) rep(paste("Group ", y), niter)))))
            plot_df[,1:2] <- as.data.frame(cbind(rep(1:niter, ngr + 1), nobs_v))
            plot_df <- rbind(plot_df[-c(1:niter),], plot_df[c(1:niter),])
            colnames(plot_df) <- c("V1", "V2", "V3")
            ggplot2::ggplot(plot_df, aes(x = V1, y = V2, color = V3)) + 
              theme_bw() +
              theme(axis.ticks = element_blank(),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank()) + 
              geom_line() +
              facet_wrap(~ factor(V3), ncol = 1) +  guides(fill=FALSE, color=FALSE)
            
            # geom_line(aes(x = V1, y = V2), size=.5, col = col) + 
            # geom_line(aes(x = V1, y = V3), size=.5, linetype="dashed", col = col) +
            # geom_line(aes(x = V1, y = V4), size=.5, linetype="dashed", col = col) 
          })

# 
# setMethod(f = "acfplot",
#           signature(x = "modCondDep"),
#           definition = function(x, maxlag = 30, alpha = 0.95){
#             stopifnot(class(x) == "modCond")
#             
#             acf_temp <- acf(x = apply(x@clust, 1, function(y) length(unique(y))), lag.max = maxlag, plot = FALSE)
#             conf.lims <- c(-1,1)*qnorm((1 + alpha)/2)/sqrt(acf_temp$n.used)
#             plot_db <- data.frame(acf_temp$lag, acf_temp$acf)
#             ggplot2::ggplot(plot_db, aes(x=acf_temp.lag, y = acf_temp.acf)) + scale_x_continuous(breaks=seq(0,maxlag + 1,round(maxlag / 4))) +
#               theme_bw() +
#               geom_hline(yintercept=conf.lims, lty=2, col='blue') +
#               labs(y="Autocorrelations", x="Lag", title= "") +
#               geom_segment(aes(xend=acf_temp.lag, yend=0)) + geom_point()
#             
#           })