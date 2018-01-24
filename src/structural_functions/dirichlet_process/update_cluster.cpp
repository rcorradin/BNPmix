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

#include <distributions/gaussian.hpp>
#include <distributions/rintnunif.hpp>
#include <distributions/tstudent.hpp>
#include <distributions/wishart.hpp>

// [[Rcpp::depends("RcppArmadillo")]]

/* 
  Update cluster, based on marginal approach of
  Dirichlet process mixture modelling. 

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
    - theta    double, precision parameter of the Dirichlet process
    - napprox: int, number of approximation for the probability of new cluster
  
  Void function.
*/

//[[Rcpp::export]]
void update_cluster_cpp(arma::mat data, 
                        arma::cube& Lambda, 
                        arma::mat& mu, 
                        arma::vec& clust, 
                        arma::vec &useful, 
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
    
    int k = (int) arma::sum(useful);
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
    
    clust[i] = rintnunif(prob, k+1);
    if(clust[i] == k){
      useful[k] = 1;
      
      arma::vec cdata = arma::trans(data.row(i)) - m0;
      arma::mat tempMat = arma::inv(sigma + cdata * arma::trans(cdata));
      Lambda.slice(k) = arma::inv(rWishartMat(nu0 + 1, tempMat));
      
      arma::mat Bn = arma::inv(arma::inv(B0) + arma::inv(Lambda.slice(k)));
      arma::vec mn = Bn * (arma::inv(B0) * m0 + arma::inv(Lambda.slice(k)) * arma::trans(data.row(i)));
      mu.row(k) = rmvnormMat(1, mn, Bn);
    }
  }
}