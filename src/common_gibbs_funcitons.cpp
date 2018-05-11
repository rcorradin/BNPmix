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
  double n = (double) clust.n_elem;

  arma::vec temp_out(grid_l);
  temp_out.fill(0);

  // for each different component
  for(int j = 0; j < k; j++){

    // evaluated the density weigthed by the frequence of the component
    temp_out += arma::sum(clust == j) * dmvnrm_ar_mat(grid, mu.row(j), Lambda.slice(j));
  }

  return temp_out / n;
}
