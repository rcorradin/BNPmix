#' @name DDPdensity
#' @export DDPdensity
#'
#' @title MCMC for GM-dependent Dirichlet process mixtures of Gaussians
#' @description The \code{DDPdensity} function generates posterior density samples for a univariate Griffiths-Milne dependent Dirichlet process mixture model with Gaussian
#' kernel, for partially exchangeable data. The function implements the importance conditional sampler method.
#'
#' @param y a vector or matrix giving the data based on which densities are to be estimated;
#' @param group vector of length \code{length(y)} containing the group labels (integers)
#' for the elements of \code{y};
#' @param output a list of arguments for generating posterior output. It contains:
#'
#' \itemize{
#' \item \code{grid}, a grid of points at which to evaluate the estimated posterior mean densities (common for all the groups).
#' \item \code{out_type}, if \code{out_type = "FULL"}, return the estimated partitions and the realizations of the posterior density for each iterations. If \code{out_type = "MEAN"}, return
#' the estimated partitions and the mean of the densities sampled at each iterations. If \code{out_type = "CLUST"}, return the estimated partitions. Default \code{out_type = "FULL"}.
# \item \code{mcmc_dens}, if equal to \code{TRUE}, the function returns a total of \code{niter}-\code{nburn} realizations of the posterior
# densities, that is one per stored iteration, evaluated at
# \code{grid} (default is \code{TRUE}). See \code{value} for details.
#' }
#'
#' @param mcmc list of MCMC arguments:
#' \itemize{
#' \item \code{niter}  (mandatory), number of iterations.
#'
#' \item \code{nburn}  (mandatory),  number of iterations to discard as burn-in.
#'
#'   \item \code{nupd}, argument controlling the number of iterations to be displayed on screen: the function reports
#'   on standard output every time \code{nupd} new iterations have been carried out (default is \code{niter/10}).
#'
#'   \item \code{print_message}, control option. If equal to \code{TRUE}, the status is printed
#'   to standard output every \code{nupd} iterations (default is \code{TRUE}).
#'
#'   \item \code{m_imp}, number of generated values for the importance sampling step of the
#'   importance conditional sampler (default is 10). See \code{details}.
#'
#' }
#'
#' @param prior a list giving the prior information, which contains:
#'
#' \itemize{
#' \item \code{strength}, the strength parameter, or total mass, of the marginal Dirichlet processes (default 1);
#' \item \code{m0}, mean of the normal base measure on the location parameter (default is the sample mean of the data);
#' \item \code{k0}, scale factor appearing in the normal base measure on the location parameter (default 1);
#' \item \code{a0}, shape parameter of the inverse gamma base measure on the scale parameter (default 2);
#' \item \code{b0}, scale parameter of the inverse gamma base measure on the scale parameter (default is the sample variance of the data);
#' \item \code{wei}, parameter controlling the strength of dependence across Dirichlet processes (default 1/2).
#' }
#'
#' @return A \code{BNPdensity} class object containing the estimated densities for each iteration,
#' the allocations for each iteration; the grid used to evaluate the densities (for each group); the
#' densities sampled from the posterior distribution (for each group); the groups; the weights of the processes.
#' The function returns also informations regarding the estimation: the number of iterations, the number
#' of burn-in iterations and the execution time.
#'
#' @details
#'
#' This function fits a Griffiths-Milne dependent Dirichlet process (GM-DDP) mixture
#' for density estimation for partially exchangeable data (Lijoi et al., 2014).
#' For each observation the \code{group} variable allows the observations to be gathered
#' into \eqn{L}=\code{length(unique(group))} distinct groups.
#' The model assumes exchangeability within each group, with observations in the \eqn{l}th group marginally
#' modelled by a location-scale Dirichlet process mixtures, i.e.
#' \deqn{\tilde f_l(y) = \int \phi(y; \mu, \sigma^2) \tilde p_l (d \mu, d \sigma^2)}
#' where each \eqn{\tilde p_l} is a Dirichlet process with total mass \code{strength} and base measure \eqn{P_0}.
#' The vector \eqn{\tilde p = (\tilde p_1,\ldots,\tilde p_L)} is assumed to be jointly distributed as a vector of
#' GM-DDP(\code{strength}, \code{wei};  \eqn{P_0}), where \code{strength} and
#' \eqn{P_0} are the total mass parameter and the base measure of each \eqn{\tilde p_l}, and \code{wei}
#' controls the dependence across the components of
#' \eqn{\tilde p}. Admissible values for \code{wei} are in \eqn{(0,1)}, with the two extremes of the range
#' corresponding to full exchangeability (\code{wei}\eqn{\rightarrow 0})
#'  and  independence across groups (\code{wei}\eqn{\rightarrow 1}).
#'
#'  \eqn{P_0} is a normal-inverse gamma base measure, i.e.
#'   \deqn{P_0(d\mu,d\sigma^2) = N(d \mu; m_0, \sigma^2 / k_0) \times IGa(d \sigma^2; a_0, b_0).}{%
#'       P_0 (d\mu,d\sigma^2) = N(d \mu; m0, \sigma^2 / k0)  IGa(d \sigma^2; a0, b0).}
#'
#' Posterior sampling is obtained by implementing the importance conditional sampler (Canale et al., 2019).
#'
#' @examples
#' data_toy <- c(rnorm(50, -4, 1), rnorm(100, 0, 1), rnorm(50, 4, 1))
#' group_toy <- c(rep(1,100), rep(2,100))
#' grid <- seq(-7, 7, length.out = 50)
#' est_model <- DDPdensity(y = data_toy, group = group_toy,
#' mcmc = list(niter = 200, nburn = 100, napprox_unif = 50),
#' output = list(grid = grid))
#' summary(est_model)
#' plot(est_model)
#'
#' @references
#' Lijoi, A., Nipoti, B., and Pruenster, I. (2014). Bayesian inference with
#' dependent normalized completely random measures. Bernoulli 20, 1260–1291.
#'
#' Canale, A., Corradin, R., & Nipoti, B. (2019). Importance conditional sampling for
#'  Bayesian nonparametric mixtures. arXiv preprint arXiv:1906.08147.


DDPdensity <- function(y,
                       group,
                       mcmc = list(),
                       prior = list(),
                       output = list()){

  if(!is.vector(y)) stop("Wrong data dimension")
  if(!is.vector(group) & !is.factor(group)) stop("Wrong group dimension")
  if(is.null(mcmc$niter)) stop("Missing number of iterations")
  if(is.null(mcmc$nburn)) mcmc$nburn = 0

  if(!is.list(mcmc)) stop("mcmc must be a list")
  if(!is.list(prior)) stop("prior must be a list")
  if(!is.list(output)) stop("output must be a list")

  if(!is.null(mcmc$niter) && (!is.numeric(mcmc$niter) | (mcmc$niter<1))) stop("mcmc$iter must be a positive integer")
  if(!is.null(mcmc$nburn) && (!is.numeric(mcmc$nburn) | (mcmc$nburn<1)) & (mcmc$nburn>mcmc$niter)) stop("mcmc$nburn must be a positive integer less than niter")
  if(!is.null(mcmc$nupd) && (!is.numeric(mcmc$nupd)  | (mcmc$nupd<1))) stop("mcmc$nupd must be a positive integer")
  if(!is.null(mcmc$m_imp) && (!is.numeric(mcmc$m_imp) | (mcmc$m_imp<1))) stop("mcmc$m_imp must be a positive integer")
  if(!is.null(mcmc$print_message) & (!is.logical(mcmc$print_message))) stop("mcmc$print_message must be a logical value")
  if(!is.null(mcmc$napprox_unif) && ((!is.numeric(mcmc$napprox_unif) | (mcmc$napprox_unif<1)))) stop("mcmc$napprox_unif must be a positive integer")


  if(!is.null(prior$m0) & !is.numeric(prior$m0)) stop("prior$m0 must be a numerical value")
  if(!is.null(prior$k0) && (!is.numeric(prior$k0) | (prior$k0<=0))) stop("prior$k0 must be a numerical positive value")
  if(!is.null(prior$a0) && (!is.numeric(prior$a0) | (prior$a0<=0))) stop("prior$a0 must be a numerical positive value")
  if(!is.null(prior$b0) && (!is.numeric(prior$b0) | (prior$b0<=0))) stop("prior$b0 must be a numerical positive value")
  if(!is.null(prior$strength) & !is.numeric(prior$strength)) stop("prior$strength must be a numerical value")
  if(!is.null(prior$wei) && (!is.numeric(prior$wei) | (prior$wei<0) | (prior$wei>1))) stop("prior$wei must be a numerical value between 0 and 1")
  if(!is.null(output$grid) & (!is.vector(output$grid))) stop("output$grid must be a vector")

  # if mcmc misses some parts, add default
  niter = mcmc$niter
  nburn = mcmc$nburn
  nupd = ifelse(is.null(mcmc$nupd), round(niter / 10), mcmc$nupd)
  print_message = ifelse(is.null(mcmc$print_message), TRUE, mcmc$print_message)
  napprox_unif = ifelse(is.null(mcmc$napprox_unif), 100, mcmc$napprox_unif)
  m_imp = ifelse(is.null(mcmc$m_imp), 10, mcmc$m_imp)

  # output
  output$out_type = ifelse(is.null(output$out_type), "FULL", output$out_type)
  if(output$out_type == "FULL"){
    mean_dens = FALSE
    mcmc_dens = TRUE
    if(is.null(output$grid)){
      grid_use = seq(from = min(y) - 0.1 * diff(range(y)), to = max(y) + 0.1 * diff(range(y)), length.out = 30)
    } else {
      if(length(dim(output$grid)) > 1) stop("Wrong grid dimension")
      grid_use <- as.vector(output$grid)
    }
  } else if (output$out_type == "MEAN"){
    mean_dens = TRUE
    mcmc_dens = TRUE
    if(is.null(output$grid)){
      grid_use = seq(from = min(y) - 0.1 * diff(range(y)), to = max(y) + 0.1 * diff(range(y)), length.out = 30)
    } else {
      if(length(dim(output$grid)) > 1) stop("Wrong grid dimension")
      grid_use <- as.vector(output$grid)
    }
  } else if (output$out_type == "CLUST"){
    mean_dens = FALSE
    mcmc_dens = FALSE
    grid_use = seq(from = min(y) - 0.1 * diff(range(y)), to = max(y) + 0.1 * diff(range(y)), length.out = 30)
  }

  group <- as.numeric(as.factor(group))
  ngr   <- length(unique(group))

  if(is.null(prior$m0)){ m0 = mean(y) } else { m0 = prior$m0 }
  if(is.null(prior$k0)){ k0 = 1 } else { k0 = prior$k0 }
  if(is.null(prior$a0)){ a0 = 2 } else { a0 = prior$a0 }
  if(is.null(prior$b0)){ b0 = var(y) } else { b0 = prior$b0 }

  if(is.null(prior$strength)){ strength = 1 } else { strength = prior$strength }
  if(is.null(prior$wei)){ wei = 0.5 } else { wei = prior$wei }


  est_model <- cDDP(y,
                    group,
                    ngr,
                    grid_use,
                    niter,
                    nburn,
                    m0,
                    k0,
                    a0,
                    b0,
                    strength,
                    wei,
                    m_imp,
                    napprox_unif,
                    nupd,
                    mcmc_dens,
                    print_message,
                    mean_dens)

  if(isTRUE(mcmc_dens)){
    if(mean_dens == TRUE){
      result <- BNPdens(density = est_model$dens[,,1],
                        grideval = grid_use,
                        clust = (est_model$clust + 1),
                        group_log = est_model$group_log,
                        niter = niter,
                        nburn = nburn,
                        tot_time = est_model$time,
                        group = group,
                        wvals = est_model$wvals,
                        dep = TRUE)
    } else {
      result <- BNPdens(density = est_model$dens,
                        grideval = grid_use,
                        clust = (est_model$clust + 1),
                        group_log = est_model$group_log,
                        niter = niter,
                        nburn = nburn,
                        tot_time = est_model$time,
                        group = group,
                        wvals = est_model$wvals,
                        dep = TRUE)
    }

  }else{
    result <- BNPdens(clust = (est_model$clust + 1),
                      group_log = est_model$group_log,
                      niter = niter,
                      nburn = nburn,
                      tot_time = est_model$time,
                      group = group,
                      wvals = est_model$wvals,
                      dep = TRUE)
  }

  return(result)
}

