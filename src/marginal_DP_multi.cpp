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

#include "para_cleanser.hpp"
#include "update_cluster.hpp"
#include "update_parameters.hpp"
#include "update_hyperparameters.hpp"
#include "update_distribution.hpp"

// [[Rcpp::depends("RcppArmadillo")]]

/*
  Update hyperparameters, according to the posterior distirbutions

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
Rcpp::List marginal_DP_multi(int nsim,
                    int nburn,
                    int napprox,
                    int nparam,
                    int d,
                    int grid_l,
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
                    arma::mat M1,
                    int s1,
                    arma::mat S1,
                    double t1,
                    double t2,
                    int nupd,
                    int plim,
                    bool FIX) {
  int n = data.n_rows;

  // initialize results
  arma::mat result_clust(nsim - nburn, n);
  arma::vec distribution(grid_l);
  arma::vec result_theta(nsim - nburn);

  // initialize cluster and parameters
  arma::mat mu(nparam, d);
  arma::cube Lambda(d,d,nparam);
  arma::vec clust(n);
  arma::vec useful(nparam);

  // fill the vectors
  distribution.fill(0);
  mu.fill(0);
  Lambda.fill(0);
  useful.fill(0);

  clust = conf_start;

  // initialize the useful parameters
  arma::vec uniqsub = unique(clust);
  int n_uniq = uniqsub.n_elem;

  for(int i = 0; i < n_uniq; i++) {
    mu.row(i) = arma::trans(mu_start);
    Lambda.slice(i) = Lambda_start;
    useful(i) = 1;
  }

  para_cleanser(Lambda = Lambda,
                mu = mu,
                clust = clust,
                useful = useful);

  int start_s = clock();
  int current_s;
  // strarting loop
  for(int sim = 0; sim < nsim; sim++){

    // update cluster allocation
    update_cluster_cpp(data = data,
                       Lambda = Lambda,
                       mu = mu,
                       clust = clust,
                       useful = useful,
                       m0 = m0,
                       B0 = B0,
                       nu0 = nu0,
                       sigma = sigma,
                       theta = theta,
                       napprox = napprox);

    // clean parameter objects
    para_cleanser(Lambda = Lambda,
                  mu = mu,
                  clust = clust,
                  useful = useful);

    // update parameters objects
    update_parameters(data = data,
                      Lambda = Lambda,
                      mu = mu,
                      clust = clust,
                      useful = useful,
                      m0 = m0,
                      B0 = B0,
                      nu0 = nu0,
                      sigma = sigma);

    // update hyperparameters
    update_hyperparameters(n = n,
                           theta = theta,
                           Lambda = Lambda,
                           mu = mu,
                           clust = clust,
                           useful = useful,
                           m0 = m0,
                           B0 = B0,
                           nu0 = nu0,
                           sigma = sigma,
                           b1 = b1,
                           B1 = B1,
                           m1 = m1,
                           M1 = M1,
                           s1 = s1,
                           S1 = S1,
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
                                                useful = useful,
                                                clust = clust,
                                                theta = theta,
                                                n = n);
    }
    if((sim + 1) % nupd == 0){
      current_s = clock();
      Rcpp::Rcout << "Completed:\t" << (sim + 1) << "/" << nsim << " - in " <<
        double(current_s-start_s)/CLOCKS_PER_SEC << " sec\n";
    }
    Rcpp::checkUserInterrupt();
  }
  distribution /= (nsim - nburn);
  Rcpp::Rcout << "Heya! Your estimation is done Captain, ARGH!\n";

  Rcpp::List resu;
  resu["dist"]  = distribution;
  resu["clust"] = result_clust;
  resu["theta"] = result_theta;
  return resu;
}
