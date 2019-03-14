#' modCond class constructor
#'
#' define the modCond class, could be for univariate or multivariate estimated models
#'
#' @param density a matrix containing the density value corresponding to the grid points
#' @param grideval a set of values to evaluate the density
#' @param clust a matrix containing the corresponding clusters for each observation, for each iteration
#' @param mean values for the location parameters
#' @param sigma2 values for the scale parameters
#' @param probs values for the mixing measure
#' @param niter number of iterations
#' @param nburn number of burn-in iterations
#' @param nnew number of new clusters sampled at each iteration
#' @param tot_time execution time
#' @param univariate if the model is univariate (TRUE/FALSE)
#' @param dep if the model is dependent Dirichlet process (TRUE/FALSE)
#' @param group_log group allocation for dependent Dirichlet process
#' @param nclust number of cluster for dependent Dirichlet process
#' @param group true group allocation for dependent Dirichlet process
#' @param wvals values of the processes weights
#'
#' @export

modCond <- function(
  density = NULL,
  grideval = NULL,
  clust = NULL,
  mean = NULL,
  sigma2 = NULL,
  probs = NULL,
  niter = NULL,
  nburn = NULL,
  nnew = NULL,
  tot_time = NULL,
  univariate = TRUE,
  dep = FALSE,
  group_log = NULL,
  nclust = NULL,
  group = NULL,
  wvals = NULL
){
  value <- list(density = density,
                grideval = grideval,
                clust = clust,
                mean = mean,
                sigma2 = sigma2,
                probs = probs,
                niter = niter,
                nburn = nburn,
                nnew = nnew,
                tot_time = tot_time,
                univariate = univariate,
                dep = dep,
                group_log = group_log,
                nclust = nclust,
                group = group,
                wvals = wvals)
  attr(value, "class") <- "modCond"
  value
}

#' modCond summary method
#'
#' @param object an object of class modCond
#' @param ... additional arguments to be passed
#'
#' @rdname summary
#' @export

summary.modCond <- function(object, ...) {
  if(!object$dep){
    if(object$univariate){
      cat("condMCMC function call:\n",
          object$nburn, "\tburn-in iterations\n",
          object$niter, "\titerations \n",
          "Global estimation time:", round(object$tot_time, digits = 2), "seconds")
    } else {
      cat("condMCMCmv function call:\n",
          object$nburn, "\tburn-in iterations\n",
          object$niter, "\titerations \n",
          "Global estimation time:", round(object$tot_time, digits = 2), "seconds")
    }
  } else {
    cat("condDDP function call:\n",
        length(table(object$group)), "\tdifferent groups\n",
        object$nburn, "\tburn-in iterations\n",
        object$niter, "\titerations \n",
        "Global estimation time:", round(object$tot_time, digits = 2), "seconds")
  }
}

#' modCond summary method
#'
#' @param x an object of class modCond
#' @param dimension if multivariate, the two dimensions for the plot (if they are equal, a marginal plot is performed)
#' @param col the color of the lines
#' @param ... additional arguments to be passed
#'
#' @rdname plot
#' @export

plot.modCond <- function(x, dimension = c(1,2), col = "#0037c4", ...) {
  if(!x$dep){
    if(x$univariate){
      with(x,{
        if(length(dim(x$density)) == 2){
          plot_df <- as.data.frame(cbind(x$grideval, colMeans(x$density)))
        } else {
          plot_df <- as.data.frame(cbind(x$grideval, x$density))
        }

        names(plot_df) = c("V1", "V2")
        ggplot2::ggplot(plot_df, mapping = ggplot2::aes(x = V1, y = V2)) +
          ggplot2::theme_bw() +
          ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                         axis.title.x = ggplot2::element_blank(),
                         axis.title.y = ggplot2::element_blank()) +
          ggplot2::geom_line(mapping = ggplot2::aes(x = V1, y = V2), size= 1, color = col)
      })
    } else {
      with(x,{
        if(dim(x$density)[2] > 1){
          plot_df <- as.data.frame(cbind(x$grideval, colMeans(x$density)))
        } else {
          plot_df <- as.data.frame(cbind(x$grideval, x$density))
        }

        names(plot_df) = c(paste("GR", 1:ncol(x$grideval), sep = ''), "V1")

        if(dimension[1] == dimension[2]){
          plot_df_use <- aggregate(plot_df, by = list(plot_df[[dimension[1]]]), FUN = sum)
          ggplot2::ggplot(data = plot_df_use, mapping = ggplot2::aes(x = Group.1, y = V1)) +
            ggplot2::theme_bw() +
            ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                           axis.title.x = ggplot2::element_blank(),
                           axis.title.y = ggplot2::element_blank()) +
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
      })
    }
  } else {
    with(x,{
      ngr <- length(unique(x$group))
      plot_df <- as.data.frame(cbind(rep(x$grideval, ngr), as.vector(apply(x$density, c(1,2), mean)),
                                     as.vector(sapply(1:ngr, function(y) rep(paste("Group ", y), length(x$grideval)))) ))
      plot_df[,1:2] <- as.data.frame(cbind(rep(x$grideval, ngr), as.vector(apply(x$density, c(1,2), mean))))

      ggplot2::ggplot(plot_df, mapping = ggplot2::aes(x = V1, y = V2, color = V3)) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                       axis.title.x = ggplot2::element_blank(),
                       axis.title.y = ggplot2::element_blank()) +
        ggplot2::geom_line() +
        ggplot2::facet_wrap(~ factor(V3), ncol = 1) +
        ggplot2::guides(fill=FALSE, color=FALSE)
    })
  }
}
