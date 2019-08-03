#' BNPdens class constructor
#'
#' A constructor for the \code{BNPdens} class. The class \code{BNPdens} is a named list containing
#' the output of the posterior estimation of different Bayesian nonparametric mixture models via
#' different MCMC strategies implemented in \code{PYdensity},  \code{DDPdensity}, and  \code{PYregression}.
#'
#' @param density a matrix containing the density value corresponding to the grid points;
#' @param data a dataset;
#' @param grideval a set of values where to evaluate the density;
#' @param grid_x regression grid, independent variable;
#' @param grid_y regression grid, dependent variable;
#' @param clust a (\code{niter - nburn}) \eqn{\times}{x} \code{nrow(data)} matrix containing
#' the cluster labels for each observation (cols) and MCMC iteration (rows);
#' @param mean values for the location parameters;
#' @param beta coefficients for regression model (only for \code{PYregression});
#' @param sigma2 values of the scale parameters;
#' @param probs values for the mixture weights;
#' @param niter number of MCMC iterations;
#' @param nburn number of MCMC iterations to discard as burn-in;
#' @param tot_time total execution time;
#' @param univariate logical. It is \code{TRUE} if the model is univariate;
#' @param regression logical. It is \code{TRUE}  for the output of  \code{PYregression};
#' @param dep  logical. It is \code{TRUE}  for the output of  \code{DDPdensity};
#' @param group_log group allocation for each iteration (only for \code{DDPdensity});
#' @param group allocation of strata of the observations (only for \code{DDPdensity});
#' @param wvals values of the processes weights (only for \code{DDPdensity}).
#'
#' @examples
#' data_toy <- c(rnorm(100, -3, 1), rnorm(100, 3, 1))
#' grid <- seq(-7, 7, length.out = 50)
#' est_model <- PYdensity(y = data_toy, mcmc = list(niter = 100,
#'                       nburn = 10, nupd = 100), output = list(grid = grid))
#' str(est_model)
#' class(est_model)
#' @export
#'

BNPdens <- function(
  density = NULL,
  data = NULL,
  grideval = NULL,
  grid_x = NULL,
  grid_y = NULL,
  clust = NULL,
  mean = NULL,
  beta = NULL,
  sigma2 = NULL,
  probs = NULL,
  niter = NULL,
  nburn = NULL,
  tot_time = NULL,
  univariate = TRUE,
  regression = FALSE,
  dep = FALSE,
  group_log = NULL,
  group = NULL,
  wvals = NULL
){
  value <- list(density = density,
                data = data,
                grideval = grideval,
                grid_x = grid_x,
                grid_y = grid_y,
                clust = clust,
                mean = mean,
                beta = beta,
                sigma2 = sigma2,
                probs = probs,
                niter = niter,
                nburn = nburn,
                tot_time = tot_time,
                univariate = univariate,
                regression = regression,
                dep = dep,
                group_log = group_log,
                group = group,
                wvals = wvals)
  attr(value, "class") <- "BNPdens"
  value
}

#' BNPpart class constructor
#'
#' A constructor for the \code{BNPpart} class. The class \code{BNPpart} is a named list containing
#' the output of the partition estimation of different methods.
#'
#' @param partitions a matrix, each row is a visited partition;
#' @param scores a vector, each value is the score of a visited partition;
#' @param psm a matrix, posterior similarity matrix.
#'
#' @examples
#' data_toy <- c(rnorm(100, -3, 1), rnorm(100, 3, 1))
#' grid <- seq(-7, 7, length.out = 50)
#' est_model <- PYdensity(y = data_toy, mcmc = list(niter = 100,
#'                       nburn = 10, nupd = 100), output = list(grid = grid))
#' part <- partition(est_model)
#' class(part)
#'
#' @export
#'

BNPpart <- function(
  partitions = NULL,
  scores = NULL,
  psm = NULL
){
  value <- list(partitions = partitions,
                scores = scores,
                psm = psm)
  attr(value, "class") <- "BNPpart"
  value
}



# -----------------------------------------------------------------------
# SUMMARY
# -----------------------------------------------------------------------

#' BNPdens summary method
#'
#' @param object an object of class \code{BNPdens};
#' @param ... additional arguments to be passed.
#'
#' @rdname summary
#'
#' @examples
#' data_toy <- c(rnorm(100, -3, 1), rnorm(100, 3, 1))
#' grid <- seq(-7, 7, length.out = 50)
#' est_model <- PYdensity(y = data_toy, mcmc = list(niter = 100,
#'                       nburn = 10, napprox = 10), output = list(grid = grid))
#' class(est_model)
#' summary(est_model)
#'
#' @export

summary.BNPdens <- function(object, ...) {
  if(!object$dep){
    if(object$univariate){
      cat("PYdensity function call:\n",
          object$nburn, "\tburn-in iterations\n",
          object$niter, "\titerations \n",
          "Global estimation time:", round(object$tot_time, digits = 2), "seconds")
    } else {
      cat("PYdensity function call:\n",
          object$nburn, "\tburn-in iterations\n",
          object$niter, "\titerations \n",
          "Global estimation time:", round(object$tot_time, digits = 2), "seconds")
    }
  } else {
    cat("DDPdensity function call:\n",
        length(table(object$group)), "\tdifferent groups\n",
        object$nburn, "\tburn-in iterations\n",
        object$niter, "\titerations \n",
        "Global estimation time:", round(object$tot_time, digits = 2), "seconds")
  }
}

# -----------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------
#' Density plot for \code{BNPdens} class
#' @description The function \code{plot.BNPdens} returns different plots for \code{BNPdens}
#' objects produced by different functions. See details.
#'
#' @details If \code{x} is a \code{BNPdens} object generated by \code{PYdensity}, the function returns
#' the posterior mean univarite or bivariate density plot.
#' If \code{x} is a \code{BNPdens} object generated by \code{PYregression}, the function returns
#' the scatterplot of the data colored as in the estimated partition.
#' If \code{x} is a \code{BNPdens} object generated by \code{DDPdensity}, the function returns
#' a wrapped plot with one density for each group.
#' The plots can be tuned in several ways: for univariate densities if \code{show_hist = TRUE},
#' then the plot shows also the histogram of the data; if \code{show_points = TRUE},
#' the plot shows also the observed points along the
#' x-axis; if \code{show_points = TRUE} and \code{show_clust = TRUE}, the points are colored
#' according to the partition estimated with the \code{partition} function.
#' For multivariate densities: if  \code{show_points = TRUE},
#' the plot shows also the scatterplot of the data;
#' if \code{show_points = TRUE} and  \code{show_clust = TRUE},
#' the points are colored according to the estimated partition.
#'
#' @param x an object of class \code{BNPdens};
#' @param dimension if \code{x} has been fitted to multivariate data,
#' \code{dimensions} specifies the two dimensions for the bivariate contour plot
#' (if they are equal, a marginal univarite plot is returned);
#' @param col the color of the lines;
#' @param show_points if \code{TRUE}, the function plots also the observations, default \code{FALSE};
#' @param show_hist if \code{TRUE}, and the model is univariate, the function plots also the histogram of the data, default  \code{FALSE};
#' @param show_clust if \code{TRUE}, the function plots also the points colored with respect to the estimated partition, default  \code{FALSE};
#' @param bin_size if \code{show_hist = TRUE}, it correponds to the size of each bin,
#' default \code{range(data) / 30};
#' @param wrap_dim a two argument vector, if the model is a \code{DDPdensity} estimation result,
#' it correponds to the number of rows and columns in the plot. Default \code{c(ngroup, 1)};
#' @param xlab label of the horizontal axis;
#' @param ylab label of the vertical axis;
#' @param band if \code{TRUE}, for the univariate case and the DDP case, the plot method returns also the credible bands;
#' @param conf_level two values vector, lower and upper levels for the confidence band (default  \code{c(0.025, 0.975)});
#' @param ... additional arguments to be passed.
#'
#' @rdname plot
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' # PYdensity example
#' data_toy <- c(rnorm(100, -3, 1), rnorm(100, 3, 1))
#' grid <- seq(-7, 7, length.out = 50)
#' est_model <- PYdensity(y = data_toy,
#'  mcmc = list(niter = 200, nburn = 100, nupd = 100),
#'  output = list(grid = grid))
#' class(est_model)
#' plot(est_model)
#'
#' # PYregression example
#' x_toy <- c(rnorm(100, 3, 1), rnorm(100, 3, 1))
#' y_toy <- c(x_toy[1:100] * 2 + 1, x_toy[101:200] * 6 + 1) + rnorm(200, 0, 1)
#' grid_x <- c(0, 1, 2, 3, 4, 5)
#' grid_y <- seq(0, 35, length.out = 50)
#' est_model <- PYregression(y = y_toy, x = x_toy,
#' mcmc = list(niter = 200, nburn = 100),
#' output = list(grid_x = grid_x, grid_y = grid_y))
#' summary(est_model)
#' plot(est_model)
#'
#' # DDPdensity example
#' data_toy <- c(rnorm(50, -4, 1), rnorm(100, 0, 1), rnorm(50, 4, 1))
#' group_toy <- c(rep(1,100), rep(2,100))
#' grid <- seq(-7, 7, length.out = 50)
#' est_model <- DDPdensity(y = data_toy, group = group_toy,
#' mcmc = list(niter = 200, nburn = 100, napprox_unif = 50),
#' output = list(grid = grid))
#' summary(est_model)
#' plot(est_model)
#'
#'
#' @export

plot.BNPdens <- function(x, dimension = c(1,2), col = "#0037c4", show_points = F,
                         show_hist = F, show_clust = F, bin_size = NULL, wrap_dim = NULL,
                         xlab = "", ylab = "", band = T, conf_level = c(0.025, 0.975), ...) {
  ggplot2::theme_set(ggplot2::theme_bw())
  if(is.null(x$density)){
    with(x, {
      xlab <- ifelse(xlab == "", "group", xlab)
      ylab <- ifelse(ylab == "", "count", ylab)
      part_temp <- partition(x)
      temp_plot <- ggplot2::ggplot(data.frame(data = as.factor(part_temp$partitions[1,]))) +
                    ggplot2::geom_bar(map = ggplot2::aes(x = data, color = data, fill = data), alpha = 0.3) +
                    ggplot2::labs(x = xlab, y = ylab) +
                    ggplot2::theme(legend.position = "none")
      temp_plot
    })
  } else {
    if(isTRUE(show_clust)){show_points <- TRUE}
    if(!x$regression){
      if(!x$dep){
        if(x$univariate){

          with(x,{
            if(is.null(bin_size)) bin_size <- (range(x$data)[2] - range(x$data)[1]) / 30
            data_plot = data.frame(V1 = x$data, V2 = partition(x)$partitions[3,])

            if(!is.null(x$density)){
              if(length(dim(x$density)) == 2){
                plot_df <- data.frame(grid = x$grideval, val = colMeans(x$density))
              } else {
                plot_df <- data.frame(grid = x$grideval, val = x$density)
              }
              if(isTRUE(band) & length(dim(x$density)) == 2){
                plot_df <- data.frame(grid = x$grideval,
                                      val  = colMeans(x$density),
                                      band_low = apply(x$density, 2, function(z) quantile(z, conf_level[1])),
                                      band_up  = apply(x$density, 2, function(z) quantile(z, conf_level[2])))
              }
            }

            temp_plot <- ggplot2::ggplot(plot_df, mapping = ggplot2::aes(x = grid, y = val))

            if(show_hist){
              temp_plot <- temp_plot + ggplot2::geom_histogram(data = data_plot, ggplot2::aes(x = V1, y = ..density..),
                                                               fill = "#EFEFEF", col = "#969696", binwidth = bin_size)
            }
            if(show_clust){
              temp_plot <- temp_plot + ggplot2::geom_point(data = data_plot, ggplot2::aes(x = V1, y =0, col = as.factor(V2)))
            }
            if(isTRUE(show_points) & !isTRUE(show_clust)){
              temp_plot <- temp_plot + ggplot2::geom_point(data = data_plot, ggplot2::aes(x = V1, y =0), color = "#646464")
            }
            if(isTRUE(band)){
              temp_plot <- temp_plot + ggplot2::geom_ribbon(data = plot_df,
                                                            mapping = ggplot2::aes(x = grid, ymin = band_low, ymax = band_up),
                                                            fill = col, alpha = 0.3)
            }

            temp_plot <- temp_plot +
              ggplot2::geom_line(mapping = ggplot2::aes(x = grid, y = val), size= 1, color = col) +
              ggplot2::theme(legend.position = "none") +
              ggplot2::labs(x = xlab, y = ylab)
            temp_plot
          })

        } else {
          with(x, {
            data_plot <- data.frame(V1 = x$data[,dimension[1]], V2 = x$data[,dimension[2]], V3 = partition(x)$partitions[3,])
            if(dim(x$density)[2] > 1){
              plot_df <- as.data.frame(cbind(x$grideval, colMeans(x$density)))
            } else {
              plot_df <- as.data.frame(cbind(x$grideval, x$density))
            }

            names(plot_df) = c(paste("GR", 1:ncol(x$grideval), sep = ''), "V1")
            if(dimension[1] == dimension[2]){
              plot_df_use <- aggregate(plot_df, by = list(plot_df[[dimension[1]]]), FUN = sum)
              temp_plot <- ggplot2::ggplot(data = plot_df_use, mapping = ggplot2::aes(x = Group.1, y = V1)) +
                ggplot2::geom_line(mapping = ggplot2::aes(x = Group.1, y = V1), size= 1, color = col) +
                ggplot2::labs(x = xlab, y = ylab)
            } else {
              plot_df_use <- aggregate(plot_df, by = list(plot_df[[dimension[1]]],plot_df[[dimension[2]]]), FUN = sum)
              temp_plot <- ggplot2::ggplot(data = plot_df_use, mapping = ggplot2::aes(x = Group.1, y = Group.2, z = V1))

              if(isTRUE(show_points) & !isTRUE(show_clust)){
                temp_plot <- temp_plot + ggplot2::geom_point(data = data_plot, ggplot2::aes(x = V1, y = V2), col = "#646464")
              }

              if(isTRUE(show_clust)){
                temp_plot <- temp_plot + ggplot2::geom_point(data = data_plot, ggplot2::aes(x = V1, y = V2, col = as.factor(V3)))
              }

              temp_plot <- temp_plot +
                ggplot2::stat_contour(data = plot_df_use, mapping = ggplot2::aes(x = Group.1, y = Group.2, z = V1), bins = 10, col = col) +
                ggplot2::theme(legend.position = "none") +
                ggplot2::labs(x = xlab, y = ylab)
              temp_plot
            }
          })
        }
      } else {
        if(is.null(wrap_dim)) wrap_dim = c(length(unique(x$group)), 1)
        with(x,{
          ngr <- length(unique(x$group))
          if(isTRUE(band)){
            plot_df <- data.frame(grid = rep(x$grideval, ngr),
                                  val  = as.vector(apply(x$density, c(1,2), mean)),
                                  group = factor(rep(levels(factor(x$group)), each = length(x$grideval)),
                                                 levels = levels(factor(x$group))),
                                  band_low = as.vector(apply(x$density, c(1,2), function(z) quantile(z, conf_level[1]))),
                                  band_up  = as.vector(apply(x$density, c(1,2), function(z) quantile(z, conf_level[2]))))
          } else {
            plot_df <- data.frame(grid = rep(x$grideval, ngr),
                                  val  = as.vector(apply(x$density, c(1,2), mean)),
                                  group = factor(rep(levels(factor(x$group)), each = length(x$grideval)),
                                                 levels = levels(factor(x$group))))
          }

          temp_plot <- ggplot2::ggplot(plot_df, mapping = ggplot2::aes(x = grid, y = val, color = col)) +
            ggplot2::labs(x = xlab, y = ylab) +
            ggplot2::facet_wrap(~ factor(group), nrow = wrap_dim[1], ncol = wrap_dim[2])

          if(isTRUE(band)){
            temp_plot <- temp_plot + ggplot2::geom_ribbon(data = plot_df,
                                 mapping = ggplot2::aes(x = grid, ymin = band_low, ymax = band_up, fill = col),
                                 colour = NA, alpha = 0.3)
            }

          temp_plot <- temp_plot +
            ggplot2::geom_line() +
            ggplot2::guides(fill=FALSE, color=FALSE)
          temp_plot
        })
      }
    } else {
      with(x,{
        plot_df <- as.data.frame(x$data)
        colnames(plot_df) = c("V1", "V2")
        part <- partition(x)$partitions[3,]

        temp_plot <- ggplot2::ggplot(plot_df, mapping = ggplot2::aes(x = V1, y = V2)) +
          ggplot2::geom_point(ggplot2::aes(x = V1, y = V2, col = as.factor(part))) +
          ggplot2::guides(fill=FALSE, color=FALSE) +
          ggplot2::labs(x = xlab, y = ylab)
        temp_plot
      })
    }
  }
}
