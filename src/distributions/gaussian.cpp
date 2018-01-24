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

/*
 * Gaussian distribution functions
 * a set of tools useful to evaluate
 * and generate from a multivariate 
 * Gaussian distribution
 */

#include "RcppArmadillo.h"
// [[Rcpp::depends("RcppArmadillo")]]

const double log2pi = std::log(2.0 * M_PI);

/* 
  Multivariate Gaussian density, 
  specified with covariance matrix.
  args: 
    - x:     vector point where evaluate the density
    - mean:  mean vector
    - sigma: covariance matrix
    - logd:  if true return the log-density, default false
  return (double) density value
*/

double dmvnrm_ar(arma::rowvec x, 
                 arma::rowvec mean,  
                 arma::mat sigma, 
                 bool logd = false) { 
  int xdim = x.n_cols;
  double out;
  
  arma::mat rooti  = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum  = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
  
  arma::vec z = rooti * arma::trans(x - mean) ;    
  out         = constants - 0.5 * arma::sum(z%z) + rootisum;     
  
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

/* 
  Multivariate Gaussian density, 
  specified with precision matrix.
  args: 
    - x:      vector point where evaluate the density
    - mean:   mean vector
    - lambda: precision matrix
    - logd:   if true return the log-density, default false
  return (double) density value
*/

double dmvnrm_prec(arma::rowvec x, 
                   arma::rowvec mean,  
                   arma::mat lambda, 
                   bool logd = false) { 
  int xdim = x.n_cols;
  double out;
  
  arma::mat rooti  = trimatu(arma::chol(lambda));
  double rootisum  = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
  
  arma::vec z = rooti * arma::trans(x - mean) ;    
  out         = constants - 0.5 * arma::sum(z%z) + rootisum;     
  
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

/* 
 Multivariate Gaussian density, more than one evaluation point,
specified with covariance matrix.
args: 
- x:     matrix where evaluate the density
- mean:  mean vector
- sigma: covariance matrix
- logd:  if true return the log-density, default false
return (vector) density value
*/

arma::vec dmvnrm_ar_mat(arma::mat x, 
                 arma::rowvec mean,  
                 arma::mat sigma, 
                 bool logd = false) { 
  
  int n = x.n_rows;
  int xdim = x.n_cols;
  arma::vec out(n);
  arma::mat rooti  = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum  = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
  
  for (int i=0; i < n; i++) {
    arma::vec z = rooti * arma::trans( x.row(i) - mean) ;    
    out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;     
  }  
  
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

/* 
  Multivariate Gaussian density, more than one evaluation point,
  specified with precision matrix.
  args: 
    - x:      matrix of points where evaluate the density
    - mean:   mean vector
    - lambda: precision matrix
    - logd:   if true return the log-density, default false
  return (vector) vector of density evaluations
*/

arma::vec dmvnrm_prec_mat(arma::mat x,  
                          arma::rowvec mean,  
                          arma::mat lambda, 
                          bool logd = false) { 
  int n = x.n_rows;
  int xdim = x.n_cols;
  arma::vec out(n);
  arma::mat rooti  = trimatu(arma::chol(lambda));
  double rootisum  = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
  
  for (int i=0; i < n; i++) {
    arma::vec z = rooti * arma::trans( x.row(i) - mean) ;    
    out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;     
  }  
  
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

/* 
  Multivariate Gaussian density, more than one evaluation point,
  passing the determinant as argument,
  specified with precision matrix.
  args: 
    - x:      matrix of points where evaluate the density
    - mean:   mean vector
    - lambda: precision matrix
    - det:    determinant of covariance matrix
  return (double) density value
*/

double dmn_prec_det(arma::rowvec x, 
                    arma::rowvec mean, 
                    arma::mat lambda, 
                    double det) {
  int dim = x.n_cols;
  double temp;
  double constant  = -(static_cast<double>(dim)/2.0) * log2pi;
  arma::rowvec xc  = (x - mean);
  temp = arma::sum(xc % (xc * lambda)); 
  return(exp(constant - 0.5 * temp + log(sqrt(det))));
}

/* 
  Multivariate Gaussian random generator, 
  specified with covariance matrix.
  NB based on "arma::randn" function.
  args: 
    - n:     number of values to generate
    - mean:  mean vector
    - sigma: covariance matrix
  return matrix of simulated values: 
    - rows: number of simulation
    - cols: number of dimension
*/

arma::mat rmvnormMat(int n, arma::vec mu, arma::mat sigma) {
  int ncols   = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}
