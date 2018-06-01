/*
 Copyright (C) 2018 Riccardo Corradin

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 3 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
*/

#include "RcppArmadillo.h"
#include "distributions.hpp"
#include "common_gibbs_functions.hpp"
// [[Rcpp::depends("RcppArmadillo")]]

/*
Update hyperparameters, according to the posterior distirbutions

args:
- n:       int, number of observations
- Lambda:  array, each slice is a precision matrix
- mu:      matrix, each row is a mean vector
- clust:   vector, each (integer) value is the cluster of the corresp. obs.
- m0:      vector, mean of the location component of the base measure
- B0:      matrix, variance of the location component of the base measure
- nu0:     double, gdl of the scale component of the base measure
- sigma    matrix, charateristic matrix of the scale component of the base measure
- b1:      int, gdl of the distribution on B0
- B1:      matrix, charateristic matrix of the distribution on B0
- m1:      vector, mean of the distribution on m0
- M1:      matrix, precision matrix of the distribution on m0
- s1:      int, gdl of the distribution on sigma
- S1:      matrix, charateristic matrix of the distribution on sigma
- t1:      double, shape parameter of the distribution on theta
- t2:      double, scale parameter of the distribution on theta

Void function.
*/

void update_hyperparameters_dep(int n,
                                double& theta,
                                arma::cube Lambda,
                                arma::mat mu,
                                arma::vec clust,
                                arma::vec& m0,
                                arma::mat& B0,
                                double nu0,
                                arma::mat sigma,
                                int b1,
                                arma::mat B1,
                                arma::vec m1,
                                double k1,
                                double t1,
                                double t2,
                                bool FIX) {
  int k = mu.n_rows;

  // update B0
  // by the posterior distribution (Inverse Wishart)
  arma::mat temp_muc = mu - arma::repmat(arma::mean(mu, 0), k, 1);
  arma::mat B1n      = arma::inv(B1 + arma::trans(temp_muc) * temp_muc + ((k1 * k) / (k1 + k)) *
                              arma::trans(arma::mean(mu, 0) - arma::trans(m0)) * (arma::mean(mu, 0) - arma::trans(m0)));
  B0                 = arma::inv(rWishartMat(b1 + k, B1n));

  // update m0
  // by the posterior distirbution (Gaussian)
  double kn     = k1 + k;
  arma::vec m1n = (k1 * m1 +  arma::trans(sum(mu, 0))) / kn;
  m0            = arma::trans(rmvnormMat(1, m1n, B0 / kn));

  // FIX = true  - no hyperprior on the mass of DP
  // FIX = false - Gamma hyperprior on the mass of DP
  if(FIX == false){

    // update the mass of DP
    // by the posterior distirbution (Gamma)
    double t_eta1 = arma::randg(1, arma::distr_param(theta + 1, 1.0))[0];
    double t_eta2 = arma::randg(1, arma::distr_param((double) n, 1.0))[0];
    double eta = t_eta1 / (t_eta1 + t_eta2);
    double pre = (t1 + k - 1) / (t1 + k - 1 + n * (t2 - log(eta)));
    double u = arma::randu(1)[0];
    int tval = (u < 1 - pre ? 1 : 0);
    theta = arma::randg(1, arma::distr_param(t1 + k - tval, 1/(t2 - log(eta))))[0];
  }
}

/*
marginal_DP_multi_indep_dep
marginal multivariate DP mixture with independent
components base measure

args:
- nsim:         int, number of iterations
- nburn:        int, number of burnin iterations
- napprox:      int, number of approximation in cluster updating step
- d:            int, number of dimension
- grid_l:       int, number of points of grid (to evaluate the density)
- data:         matrix, given data
- grid:         matrix, points to evaluate the density
- Lambda_start: matrix, initial value for the scale component of the kernel
- mu_start:     matrix, initial value for the location component of the kernel
- clust:        vector, each (integer) value is the cluster of the corresp. obs.
- theta:        double, precision parameter of the Dirichlet process
- m0:           vector, mean of the location component of the base measure
- B0:           matrix, variance of the location component of the base measure
- nu0:          double, gdl of the scale component of the base measure
- sigma         matrix, charateristic matrix of the scale component of the base measure
- b1:           int, gdl of the distribution on B0
- B1:           matrix, charateristic matrix of the distribution on B0
- m1:           vector, mean of the distribution on m0
- M1:           matrix, precision matrix of the distribution on m0
- s1:           int, gdl of the distribution on sigma
- S1:           matrix, charateristic matrix of the distribution on sigma
- t1:           double, shape parameter of the distribution on theta
- t2:           double, scale parameter of the distribution on theta

return list:
- distribution
- MCMC cluster allocation
- theta vector
*/

//[[Rcpp::export]]
Rcpp::List marginal_DP_multi_indep_dep(int nsim,
                                       int nburn,
                                       int napprox,
                                       int d,
                                       arma::mat data,
                                       arma::mat grid,
                                       arma::vec conf_start,
                                       arma::vec mu_start,
                                       arma::mat Lambda_start,
                                       double theta,
                                       arma::vec m0,
                                       arma::mat B0,
                                       double nu0,
                                       arma::mat sigma,
                                       int b1,
                                       arma::mat B1,
                                       arma::vec m1,
                                       double k1,
                                       double t1,
                                       double t2,
                                       int nupd,
                                       bool FIX) {
  int n = data.n_rows;
  int grid_l = grid.n_rows;

  // initialize results
  arma::mat result_clust(nsim - nburn, n);
  arma::vec distribution(grid_l);
  arma::vec result_theta(nsim - nburn);
  distribution.fill(0);

  // initialize cluster
  arma::vec clust(n);
  clust = conf_start;

  // initialize the useful parameters
  arma::vec uniqsub = unique(clust);
  int n_uniq = uniqsub.n_elem;

  arma::mat mu(n_uniq, d);
  arma::cube Lambda(d,d,n_uniq);
  mu.fill(0);
  Lambda.fill(0);

  for(int i = 0; i < n_uniq; i++) {
    mu.row(i) = arma::trans(mu_start);
    Lambda.slice(i) = Lambda_start;
  }

  para_cleanser(Lambda = Lambda,
                mu = mu,
                clust = clust);

  int start_s = clock();
  int current_s;
  // strarting loop
  for(int sim = 0; sim < nsim; sim++){

    // update cluster allocation
    update_cluster_cpp(data = data,
                       Lambda = Lambda,
                       mu = mu,
                       clust = clust,
                       m0 = m0,
                       B0 = B0,
                       nu0 = nu0,
                       sigma = sigma,
                       theta = theta,
                       napprox = napprox);

    // clean parameter objects
    para_cleanser(Lambda = Lambda,
                  mu = mu,
                  clust = clust);

    // update parameters objects
    update_parameters(data = data,
                      Lambda = Lambda,
                      mu = mu,
                      clust = clust,
                      m0 = m0,
                      B0 = B0,
                      nu0 = nu0,
                      sigma = sigma);

    // update hyperparameters
    update_hyperparameters_dep(n = n,
                               theta = theta,
                               Lambda = Lambda,
                               mu = mu,
                               clust = clust,
                               m0 = m0,
                               B0 = B0,
                               nu0 = nu0,
                               sigma = sigma,
                               b1 = b1,
                               B1 = B1,
                               m1 = m1,
                               k1 = k1,
                               t1 = t1,
                               t2 = t2,
                               FIX = FIX);

    // save quantities
    if(sim >= nburn){

      result_clust.row(sim - nburn) = arma::trans(clust);
      result_theta(sim - nburn) = theta;

      // update distributrion
      distribution += update_distribution(grid = grid,
                                          grid_l = grid_l,
                                          mu = mu,
                                          Lambda = Lambda,
                                          clust = clust,
                                          theta = theta);
    }
    if((sim + 1) % nupd == 0){
      current_s = clock();
      Rcpp::Rcout << "Completed:\t" << (sim + 1) << "/" << nsim << " - in " <<
        double(current_s-start_s)/CLOCKS_PER_SEC << " sec\n";
    }
    Rcpp::checkUserInterrupt();
  }
  Rcpp::Rcout << "Heya! Your estimation is done Captain, ARGH!\n";

  distribution = distribution / (nsim - nburn);

  Rcpp::List resu;
  resu["dist"]  = distribution;
  resu["clust"] = result_clust;
  resu["theta"] = result_theta;
  return resu;
}
