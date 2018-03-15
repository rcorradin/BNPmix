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
// [[Rcpp::depends("RcppArmadillo")]]
const double log2pi = std::log(2.0 * M_PI);

//  / ___| __ _ _   _ ___ ___(_) __ _ _ __
// | |  _ / _` | | | / __/ __| |/ _` | '_ \
// | |_| | (_| | |_| \__ \__ \ | (_| | | | |
//  \____|\__,_|\__,_|___/___/_|\__,_|_| |_|

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

//  _       ____  _             _            _
// | |_    / ___|| |_ _   _  __| | ___ _ __ | |_
// | __|___\___ \| __| | | |/ _` |/ _ \ '_ \| __|
// | ||_____|__) | |_| |_| | (_| |  __/ | | | |_
//  \__|   |____/ \__|\__,_|\__,_|\___|_| |_|\__|

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

// __        ___     _                _
// \ \      / (_)___| |__   __ _ _ __| |_
//  \ \ /\ / /| / __| '_ \ / _` | '__| __|
//   \ V  V / | \__ \ | | | (_| | |  | |_
//    \_/\_/  |_|___/_| |_|\__,_|_|   \__|

/*
 Generate from a Wishart distribution
 args:
 - df:     degrees of freedom
 - lambda: charateristic matrix
 return (double) matrix generated from a Wishart distribution
 */

arma::mat rWishartMat(int df,
                      arma::mat lambda) {
  int d = lambda.n_cols;
  arma::mat Y = arma::randn(df, d) * arma::chol(lambda);
  return Y.t() * Y ;
}

//  _   _                   _   _       _  __
// | \ | | ___  _ __       | | | |_ __ (_)/ _| ___  _ __ _ __ ___
// |  \| |/ _ \| '_ \ _____| | | | '_ \| | |_ / _ \| '__| '_ ` _ \
// | |\  | (_) | | | |_____| |_| | | | | |  _| (_) | |  | | | | | |
// |_| \_|\___/|_| |_|      \___/|_| |_|_|_|  \___/|_|  |_| |_| |_|
//  ____  _                   _
// |  _ \(_)___  ___ _ __ ___| |_ ___
// | | | | / __|/ __| '__/ _ \ __/ _ \
// | |_| | \__ \ (__| | |  __/ ||  __/
// |____/|_|___/\___|_|  \___|\__\___|

int rintnunif(arma::vec prob,
              int a){
  double u = arma::randu(1)[0];
  prob = arma::cumsum(prob);

  // NB now is 0 to k with loop up to a+1
  for(int k = 0; k < a + 1; k++) {
    if(u <= prob[k]) {
      return k;
    }
  }
  return 0;
}
