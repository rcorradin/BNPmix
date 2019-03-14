#' modCond class
#'
#' @slot density a matrix containing the density value corresponding to the grid points
#' @slot grideval a set of values to evaluate the density
#' @slot clust a matrix containing the corresponding clusters for each observation, for each iteration
#' @slot mean values for the location parameters
#' @slot sigma2 values for the scale parameters
#' @slot probs values for the mixing measure
#' @slot niter number of iterations
#' @slot nburn number of burn-in iterations
#' @slot nnew number of new clusters sampled at each iteration
#' @slot tot_time execution time

modCond <- setClass(Class = "modCond",
                    slots = c(density   = "matrix",
                              grideval  = "numeric",
                              clust     = "matrix",
                              mean      = "list",
                              sigma2    = "list",
                              probs     = "list",
                              niter     = "numeric",
                              nburn     = "numeric",
                              nnew      = "numeric",
                              tot_time  = "numeric"),
                    prototype = list(density  = matrix(0,0,0),
                                     grideval = numeric(),
                                     clust    = matrix(0,0,0),
                                     mean     = list(0),
                                     sigma2   = list(0),
                                     probs    = list(0),
                                     niter    = numeric(),
                                     nburn    = numeric(),
                                     nnew     = numeric(),
                                     tot_time = numeric()))

#-----------------------------------------------------------------------

#' modCondMv class
#'
#' @slot density a matrix containing the density value corresponding to the grid points
#' @slot grideval a set of values to evaluate the density
#' @slot clust a matrix containing the corresponding clusters for each observation, for each iteration
#' @slot mean values for the location parameters
#' @slot sigma2 values for the scale parameters
#' @slot probs values for the mixing measure
#' @slot niter number of iterations
#' @slot nburn number of burn-in iterations
#' @slot nnew number of new clusters sampled at each iteration
#' @slot tot_time execution time

modCondMv <- setClass(Class = "modCondMv",
                      slots = c(density   = "matrix",
                                grideval  = "matrix",
                                clust     = "matrix",
                                mean      = "list",
                                sigma2    = "list",
                                probs     = "list",
                                niter     = "numeric",
                                nburn     = "numeric",
                                nnew      = "numeric",
                                tot_time  = "numeric",
                                # TEMP
                                risk      = "matrix"),
                      prototype = list(density  = matrix(0,0,0),
                                       grideval = matrix(0,0,0),
                                       clust    = matrix(0,0,0),
                                       mean     = list(0),
                                       sigma2   = list(0),
                                       probs    = list(0),
                                       niter    = numeric(),
                                       nburn    = numeric(),
                                       nnew     = numeric(),
                                       tot_time = numeric(),
                                       # TEMP
                                       risk     = matrix(0,0,0)))

#-----------------------------------------------------------------------

#' modCondDep class
#'
#' @slot density a matrix containing the density value corresponding to the grid points
#' @slot grideval a set of values to evaluate the density
#' @slot clust a matrix containing the corresponding clusters for each observation, for each iteration
#' @slot group_log indicator variables to trace if the observation falls into its gorup or in the common one
#' @slot niter number of iterations
#' @slot nburn number of burn-in iterations
#' @slot nclust number of cluster
#' @slot tot_time execution time
#' @slot group the prior group allocation
#' @slot wvals values for the weights of the specific processes

modCondDep <- setClass(Class = "modCondDep",
                    slots = c(density   = "array",
                              grideval  = "numeric",
                              clust     = "matrix",
                              group_log      = "matrix",
                              niter     = "numeric",
                              nburn     = "numeric",
                              nclust    = "numeric",
                              tot_time  = "numeric",
                              group     = "numeric",
                              wvals     = "matrix"),
                    prototype = list(density  = array(0,0),
                                     grideval = numeric(),
                                     clust    = matrix(0,0,0),
                                     group_log     = matrix(0,0,0),
                                     niter    = numeric(),
                                     nburn    = numeric(),
                                     nclust   = numeric(),
                                     tot_time = numeric(),
                                     group    = numeric(),
                                     wvals    = matrix(0,0,0)))

#-----------------------------------------------------------------------

