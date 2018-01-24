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
 * t-Student distribution functions
 * a set of tools useful to evaluate
 * the t-Student distribution
 */

#include "RcppArmadillo.h"
// [[Rcpp::depends("RcppArmadillo")]]

/* 
  Multivariate t-Student density, 
  specified with covariance matrix.
  args: 
    - x:     vector point where evaluate the density
    - mean:  mean vector
    - sigma: covariance matrix
    - df:    degrees of freedom
    - logd:  if true return the log-density, default false
  return (double) density value
*/

double dmvt_ar(arma::rowvec x, 
               arma::rowvec mean,  
               arma::mat sigma, 
               int df, 
               bool logd = false) { 
  int d = x.n_elem;
  double out;
  
  arma::mat rooti  = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum  = arma::sum(log(rooti.diag()));
  double c = lgamma((d + df)/2.0) - lgamma(df/2.0) - (d/2.0) * std::log(M_PI * df) + rootisum;
  
  arma::vec z = rooti * arma::trans(x - mean) ;    
  out         = c - 0.5 * (df + d)  * std::log1p(arma::sum(z%z) / df);     
  
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

/* 
  Multivariate t-Student density, 
  specified with precision matrix.
  args: 
    - x:      vector point where evaluate the density
    - mean:   mean vector
    - lambda: precision matrix
    - df:     degrees of freedom
    - logd:   if true return the log-density, default false
  return (double) density value
*/

double dmvt_prec(arma::rowvec x, 
                 arma::rowvec mean,  
                 arma::mat lambda, 
                 int df, 
                 bool logd = false) { 
  int d = x.n_elem;
  double out;
  
  arma::mat rooti  = trimatu(arma::chol(lambda));
  double rootisum  = arma::sum(log(rooti.diag()));
  double c = lgamma((d + df)/2.0) - lgamma(df/2.0) - (d/2.0) * std::log(M_PI * df) + rootisum;
  
  arma::vec z = rooti * arma::trans(x - mean) ;    
  out         = c - 0.5 * (df + d)  * std::log1p(arma::sum(z%z) / df);     
  
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

/* 
  Multivariate t-Student density, more than one evaluation point, 
  specified with precision matrix.
  args: 
    - x:      matrix point where evaluate the density
    - mean:   mean vector
    - lambda: precision matrix
    - df:     degrees of freedom
    - logd:   if true return the log-density, default false
  return (double) density value
*/

arma::vec dmvt_prec_mat(arma::mat x, 
                 arma::rowvec mean,  
                 arma::mat lambda, 
                 int df, 
                 bool logd = false) { 
  int n = x.n_rows;
  int d = x.n_cols;
  arma::vec out(n);
  
  arma::mat rooti  = trimatu(arma::chol(lambda));
  double rootisum  = arma::sum(log(rooti.diag()));
  double c = lgamma((d + df)/2.0) - lgamma(df/2.0) - (d/2.0) * std::log(M_PI * df) + rootisum;

  for (int i=0; i < n; i++) {
    arma::vec z = rooti * arma::trans( x.row(i) - mean) ;    
    out(i)      = c - 0.5 * (df + d)  * std::log1p(arma::sum(z%z) / df);    
  }  
  
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}
