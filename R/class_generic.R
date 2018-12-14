#' A class
#' @export modCondMv

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
                                tot_time  = "numeric"),
                      prototype = list(density  = matrix(0,0,0),
                                       grideval = matrix(0,0,0),
                                       clust    = matrix(0,0,0),
                                       mean     = list(0),
                                       sigma2   = list(0),
                                       probs    = list(0),
                                       niter    = numeric(),
                                       nburn    = numeric(),
                                       nnew     = numeric(),
                                       tot_time = numeric()))
#-----------------------------------------------------------------------

#' A class
#' @export modCond

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

#' A class
#' @export modCondDep

modCondDep <- setClass(Class = "modCondDep",
                    slots = c(density   = "array",
                              grideval  = "numeric",
                              clust     = "matrix",
                              zeta      = "matrix",
                              niter     = "numeric",
                              nburn     = "numeric",
                              nclust    = "numeric",
                              tot_time  = "numeric",
                              group     = "numeric",
                              wvals     = "matrix",
                              vvals     = "matrix"),
                    prototype = list(density  = array(0,0),
                                     grideval = numeric(),
                                     clust    = matrix(0,0,0),
                                     zeta     = matrix(0,0,0),
                                     niter    = numeric(),
                                     nburn    = numeric(),
                                     nclust   = numeric(),
                                     tot_time = numeric(),
                                     group    = numeric(),
                                     wvals    = matrix(0,0,0),
                                     vvals    = matrix(0,0,0)))

#-----------------------------------------------------------------------

setGeneric("trace_ngr", function(x, ...) standardGeneric("trace_ngr"))
setGeneric("trace_obs", function(x, ...) standardGeneric("trace_obs"))
