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

//[[Rcpp::export]]
void update_parameters(arma::mat data, 
                       arma::cube& Lambda, 
                       arma::mat& mu, 
                       arma::vec& clust, 
                       arma::vec useful, 
                       arma::vec m0, 
                       arma::mat B0, 
                       double nu0, 
                       arma::mat sigma) {
  int k = (int) arma::sum(useful);

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