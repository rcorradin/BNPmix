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
  Update hyperparameters, according to the posterior distirbutions

  args:
    - n:       int, number of observations
    - Lambda:  array, each slice is a precision matrix
    - mu:      matrix, each row is a mean vector
    - clust:   vector, each (integer) value is the cluster of the corresp. obs.
    - useful:  vector, binary. Value 1 implies that the corresponding position is
                              used, value 0 implies that is not used
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

//[[Rcpp::export]]
double update_hyperparameters(int n, 
                              double& theta, 
                              arma::cube Lambda, 
                              arma::mat mu, 
                              arma::vec clust, 
                              arma::vec useful, 
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
                              double t2) {
  int k = (int) arma::sum(useful);
  
  arma::mat temp_mu  = mu.rows(arma::find(useful == 1));
  arma::mat temp_muc = temp_mu - arma::repmat(arma::trans(m0), k, 1);
  arma::mat B1n      = arma::inv(B1 + arma::trans(temp_muc) * temp_muc);
  B0                 = arma::inv(rWishartMat(b1 + k, B1n));
  
  // arma::mat M1n = arma::inv(M1 + k * B0);
  arma::mat M1n = arma::inv(arma::inv(M1) + k * arma::inv(B0));
  arma::vec m1n = M1n * (arma::inv(M1) * m1 + arma::inv(B0) * arma::trans(sum(temp_mu, 0)));  
  m0            = arma::trans(rmvnormMat(1, m1n, M1n));
  
  // arma::mat tMatL = arma::sum(Lambda.slices(0, k - 1), 2);
  arma::cube tinv = Lambda.slices(0, k - 1);
  for(int ind = 0; ind < k; ind++){
    tinv.slice(ind) = arma::inv(Lambda.slice(ind));
  }
  arma::mat tMatL = arma::sum(tinv, 2);
  // arma::mat tLambdas = arma::inv(S1 + tMatL);
  // sigma = arma::inv(rWishartMat(s1 + k * nu0, tLambdas));
  arma::mat tLambdas = arma::inv(arma::inv(S1) + tMatL);
  sigma = rWishartMat(s1 + k * nu0, tLambdas);
  
  // generate beta distribution as ratio of gamma distributions
  // NB the gamma distribution R::rgamma
  // is specified in function of shape and scale
  
  double eta = R::rbeta(theta + 1, n);
  double pre = (t1 + k - 1) / (t1 + k - 1 + n * (t2 - log(eta)));
  double u = arma::randu(1)[0];
  int tval = (u < 1 - pre ? 1 : 0);
  theta = R::rgamma(t1 + k - tval, 1 / (t2 - log(eta)));
}
