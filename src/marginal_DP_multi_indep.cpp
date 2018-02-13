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
// [[Rcpp::depends("RcppArmadillo")]]

/*
   Clean parameter, discard the middle not used values for the clusters and
   update the correspondent parameters.

   args:
   - Lambda:  array, each slice is a precision matrix
   - mu:      matrix, each row is a mean vector
   - clust:   vector, each (integer) value is the cluster of the corresp. obs.

   Void function.
*/

void para_cleanser(arma::cube &Lambda,
                   arma::mat &mu,
                   arma::vec &clust) {
  int k = mu.n_rows;

  // for all the used parameters
  for(int i = 0; i < k; i++){

    // if a cluster is empty
    if((int) arma::sum(clust == i) == 0){

      // find the last full cluster, then swap
      for(int j = k; j > i; j--){
        if((int) arma::sum(clust == j) != 0){
          clust( arma::find(clust == j) ).fill(i);
          mu.swap_rows(i,j);
          Lambda.slice(i).swap(Lambda.slice(j));
          break;
        }
      }
    }
  }

  // reduce dimensions
  int u_bound = 0;
  for(int i = 0; i < k; i++){
    if(arma::sum(clust == i) > 0){
      u_bound += 1;
    }
  }
  mu.resize(u_bound, mu.n_cols);
  Lambda.resize(Lambda.n_rows, Lambda.n_cols, u_bound);
}

/*
  Update cluster, based on marginal approach of
  Dirichlet process mixture modelling.

  args:
  - data:    matrix of given data
  - Lambda:  array, each slice is a precision matrix
  - mu:      matrix, each row is a mean vector
  - clust:   vector, each (integer) value is the cluster of the corresp. obs.
  - m0:      vector, mean of the location component of the base measure
  - B0:      matrix, variance of the location component of the base measure
  - nu0:     double, gdl of the scale component of the base measure
  - sigma    matrix, charateristic matrix of the scale component of the base measure
  - theta    double, precision parameter of the Dirichlet process
  - napprox: int, number of approximation for the probability of new cluster

  Void function.
*/

void update_cluster_cpp(arma::mat data,
                        arma::cube &Lambda,
                        arma::mat &mu,
                        arma::vec &clust,
                        arma::vec m0,
                        arma::mat B0,
                        double nu0,
                        arma::mat sigma,
                        double theta,
                        int napprox) {

  // initialize the number of observation

  int n = data.n_rows;

  // for each observation:
  //    using a marginal approach to generate the new
  //    cluster, given the remaining observations

  for(int i = 0; i < n; i++) {

    bool req_clean = false;
    if(arma::sum(clust == clust[i]) == 1){
      req_clean = true;
    }

    int k = mu.n_rows;
    arma::vec prob(k+1);
    prob.fill(0);
    arma::vec temp_clust = clust;
    temp_clust(i) = k+1;

    // for each cluster:
    //    compute the probability to fall inside

    for(int j = 0; j < k; j++) {
      int nj = (int) arma::sum(temp_clust == j);
      prob[j] = dmvnrm_ar(data.row(i), mu.row(j), Lambda.slice(j)) * nj;
    }

    // compute the probability to generate a new cluster

    double temp = 0.0;
    arma::vec cdata   = arma::trans(data.row(i)) - m0;
    arma::mat tempMat = arma::inv(sigma + cdata * arma::trans(cdata));
    for(int a = 0; a < napprox; a++) {
      arma::mat tWish = arma::inv(rWishartMat(nu0 + 1, tempMat));  // CONTROLLA SE LE INVERSE SONO CORRETTE
      temp += dmvnrm_ar(data.row(i), arma::trans(m0), tWish);
    }
    prob[k] = theta * (temp / napprox);
    prob = prob / arma::sum(prob);

    // generate the cluster
    // plus acceleration step:
    //    if the cluster is new, then generate a
    //    new value for the parameters

    clust[i] = rintnunif(prob, k);
    if(clust[i] == k){

      // resize the parameters objects
      // and propose the new values
      arma::vec cdata = arma::trans(data.row(i)) - m0;
      arma::mat tempMat = arma::inv(sigma + cdata * arma::trans(cdata));
      Lambda.resize(Lambda.n_rows, Lambda.n_cols, k + 1);
      Lambda.slice(k) = arma::inv(rWishartMat(nu0 + 1, tempMat));

      arma::mat Bn = arma::inv(arma::inv(B0) + arma::inv(Lambda.slice(k)));
      arma::vec mn = Bn * (arma::inv(B0) * m0 + arma::inv(Lambda.slice(k)) * arma::trans(data.row(i)));
      mu.resize(k + 1, mu.n_cols);
      mu.row(k) = rmvnormMat(1, mn, Bn);
    }

    // if required, cleans the parameters
    if(req_clean){
      para_cleanser(Lambda,
                    mu,
                    clust);
    }
  }
}

/*
  Update parameters, according to the posterior distirbutions

  args:
  - data:    matrix of given data
  - Lambda:  array, each slice is a precision matrix
  - mu:      matrix, each row is a mean vector
  - clust:   vector, each (integer) value is the cluster of the corresp. obs.
  - useful:  vector, binary. Value 1 implies that the corresponding position is
  used, value 0 implies that is not used
  - m0:      vector, mean of the location component of the base measure
  - B0:      matrix, variance of the location component of the base measure
  - nu0:     double, gdl of the scale component of the base measure
  - sigma    matrix, charateristic matrix of the scale component of the base measure

  Void function.
*/

void update_parameters(arma::mat data,
                       arma::cube &Lambda,
                       arma::mat &mu,
                       arma::vec &clust,
                       arma::vec m0,
                       arma::mat B0,
                       double nu0,
                       arma::mat sigma) {
  int k = mu.n_rows;

  // for the parameter in using
  for(int j = 0; j < k; j++){

    // update the scale of each component
    int nj = arma::sum(clust == j);
    arma::mat temp_data  = data.rows(arma::find(clust == j));
    arma::mat temp_datac = temp_data - arma::repmat(mu.row(j), nj, 1);
    arma::mat temp_mat   = arma::inv(sigma + arma::trans(temp_datac) * temp_datac);
    Lambda.slice(j)      = arma::inv(rWishartMat(nu0 + nj, temp_mat));

    // update the location of each component
    arma::mat Bn = arma::inv(arma::inv(B0) + nj * arma::inv(Lambda.slice(j)));
    arma::vec mn = Bn * (arma::inv(B0) * m0 + arma::inv(Lambda.slice(j)) * arma::trans(sum(temp_data, 0)));
    mu.row(j)    = rmvnormMat(1, mn, Bn);
  }
}

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

double update_hyperparameters(int n,
                              double& theta,
                              arma::cube Lambda,
                              arma::mat mu,
                              arma::vec clust,
                              arma::vec& m0,
                              arma::mat& B0,
                              double nu0,
                              arma::mat& sigma,
                              int b1,
                              arma::mat B1,
                              arma::vec m1,
                              arma::mat M1,
                              int s1,
                              arma::mat S1,
                              double t1,
                              double t2,
                              bool FIX) {
  int k = mu.n_rows;

  // update B0
  // by the posterior distribution (Inverse Wishart)
  arma::mat temp_muc = mu - arma::repmat(arma::trans(m0), k, 1);
  arma::mat B1n      = arma::inv(B1 + arma::trans(temp_muc) * temp_muc);
  B0                 = arma::inv(rWishartMat(b1 + k, B1n));

  // update m0
  // by the posterior distirbution (Gaussian)
  arma::mat M1n = arma::inv(arma::inv(M1) + k * arma::inv(B0));
  arma::vec m1n = M1n * (arma::inv(M1) * m1 + arma::inv(B0) * arma::trans(sum(mu, 0)));
  m0            = arma::trans(rmvnormMat(1, m1n, M1n));

  // update sigma
  // by the posterior distirbution (Wishart)
  arma::cube tinv = Lambda;
  for(int ind = 0; ind < k; ind++){
    tinv.slice(ind) = arma::inv(Lambda.slice(ind));
  }
  arma::mat tMatL = arma::sum(tinv, 2);
  arma::mat tLambdas = arma::inv(arma::inv(S1) + tMatL);
  sigma = rWishartMat(s1 + k * nu0, tLambdas);

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
   Update distribution (approximate).

   args:
   - grid:    matrix, given points to evaluate the density
   - Lambda:  array, each slice is a precision matrix
   - mu:      matrix, each row is a mean vector
   - clust:   vector, each (integer) value is the cluster of the corresp. obs.
   - useful:  vector, binary. Value 1 implies that the corresponding position is
   used, value 0 implies that is not used
   - theta    double, precision parameter of the Dirichlet process
   - n:       int, number of observation

   Void function.
*/

arma::vec update_distribution(arma::mat grid,
                              int grid_l,
                              arma::mat mu,
                              arma::cube Lambda,
                              arma::vec clust,
                              double theta){
  int k = mu.n_rows;
  int n = clust.n_elem;

  arma::vec temp_out(grid_l);
  temp_out.fill(0);

  // for each different component
  for(int j = 0; j < k; j++){

    // evaluated the density weigthed by the frequence of the component
    temp_out += arma::sum(clust == j) * dmvnrm_ar_mat(grid, mu.row(j), Lambda.slice(j));
  }

  return temp_out / n;
}

/*
  marginal_DP_multi_intep
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
Rcpp::List marginal_DP_multi_indep(int nsim,
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
                                   arma::mat M1,
                                   int s1,
                                   arma::mat S1,
                                   double t1,
                                   double t2,
                                   int nupd,
                                   bool FIX) {
  int n = data.n_rows;
  int grid_l = grid.n_rows;

  // initialize results
  arma::mat result_clust(nsim - nburn, n);
  arma::mat distribution(nsim - nburn, grid_l);
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
    update_hyperparameters(n = n,
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
      distribution.row(sim - nburn) = arma::trans(update_distribution(grid = grid,
                                                                       grid_l = grid_l,
                                                                       mu = mu,
                                                                       Lambda = Lambda,
                                                                       clust = clust,
                                                                       theta = theta));
    }
    if((sim + 1) % nupd == 0){
      current_s = clock();
      Rcpp::Rcout << "Completed:\t" << (sim + 1) << "/" << nsim << " - in " <<
        double(current_s-start_s)/CLOCKS_PER_SEC << " sec\n";
    }
    Rcpp::checkUserInterrupt();
  }
  Rcpp::Rcout << "Heya! Your estimation is done Captain, ARGH!\n";

  Rcpp::List resu;
  resu["dist"]  = distribution;
  resu["clust"] = result_clust;
  resu["theta"] = result_theta;
  return resu;
}
