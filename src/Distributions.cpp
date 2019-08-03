/*==================================================================================
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
 ==================================================================================*/

#include "RcppArmadillo.h"
// [[Rcpp::depends("RcppArmadillo")]]

/*==================================================================================
 * rintnunif - sample NON-Uniform discrete distribution
 *
 * args:
 * - weights, a vector of (eventually not normalized) probabilities
 *
 * return:
 * - an integer
 ==================================================================================*/
int rintnunif(arma::vec weights){
  double u = arma::randu();
  arma::vec probs = weights / sum(weights);
  probs = arma::cumsum(probs);

  for(arma::uword k = 0; k < probs.n_elem; k++) {
    if(u <= probs[k]) {
      return k;
    }
  }
  return -1;
}

/*==================================================================================
 * rintnunif_log - sample NON-Uniform discrete distribution - log scale
 *
 * args:
 * - lweights, a vector of log-probabilities
 *
 * return:
 * - an integer
 ==================================================================================*/
int rintnunif_log(arma::vec lweights){

  double u = arma::randu();
  arma::vec probs(lweights.n_elem);

  // for(arma::uword k = 0; k < probs.n_elem; k++) {
  //   probs(k) = 1 / sum(exp(lweights - lweights(k)));
  // }

  for(arma::uword k = 0; k < probs.n_elem; k++) {
    probs(k) = 1 / sum(exp(lweights - lweights(k)));
  }

  probs = arma::cumsum(probs);

  for(arma::uword k = 0; k < probs.n_elem; k++) {
    if(u <= probs[k]) {
      return k;
    }
  }
  return -1;
}

/*==================================================================================
 * rintnunifw - sample NON-Uniform discrete distribution (plus mass)
 *
 * args:
 * - freq, a vector of frequencies
 * - mass, mass of a process
 *
 * return:
 * - an integer
 ==================================================================================*/
int rintnunifw(arma::vec freq,
               double mass){
  arma::vec weights(freq);
  weights.resize(freq.n_elem + 1);
  weights[freq.n_elem] = mass;
  arma::vec probs = arma::cumsum(weights) / sum(weights);
  double u = arma::randu();

  for(arma::uword k = 0; k < probs.n_elem; k++) {
    if(u <= probs[k]) {
      return k;
    }
  }
  return -1;
}

/*==================================================================================
 * rdirich_mass - sample Dirichlet distribution (plus mass)
 *
 * args:
 * - freq, a vector of frequencies
 * - mass, mass of a process
 *
 * return:
 * - a vector
 ==================================================================================*/
arma::vec rdirich_mass(arma::vec freq,
                       double mass){

  arma::vec weights(freq);
  weights.resize(freq.n_elem + 1);
  weights[freq.n_elem] = mass;
  arma::vec result(weights.n_elem);

  for(arma::uword j = 0; j < weights.n_elem; j++){
    result[j] = arma::randg(1, arma::distr_param(weights[j], 1.0))[0];
  }
  return(result / sum(result));
}

/*==================================================================================
 * rdirich_mass_tot - sample Dirichlet distribution (plus mass, return total sum also)
 *
 * args:
 * - freq, a vector of frequencies
 * - mass, mass of a process
 *
 * return:
 * - a vector
 ==================================================================================*/
arma::vec rdirich_mass_tot(arma::vec freq,
                           double mass){

  arma::vec weights(freq);
  weights.resize(freq.n_elem + 1);
  weights[freq.n_elem] = mass;

  arma::vec result(weights.n_elem + 1);
  result.fill(0.0);

  for(arma::uword j = 0; j < weights.n_elem; j++){
    result[j] = arma::randg(1, arma::distr_param(weights[j], 1.0))[0];
  }

  double tot = sum(result);
  result = result / tot;
  result[result.n_elem - 1] = tot;

  return(result);
}

/*==================================================================================
 * dt_ls - evaluate univariate t-Student distribution (log scale)
 *
 * args:
 * - x, a point of the support
 * - df, degree of freedom
 * - mu, mean
 * - sigma, standard deviation
 *
 * return:
 * - double, log density
 ==================================================================================*/
double dt_ls(double x,
             double df,
             double mu,
             double sigma){
  double z = (x - mu)/sigma;
  double out = lgamma((df + 1) / 2) - log(sqrt(M_PI * df)) - log(sigma) -
    lgamma(df / 2) - (df + 1) * log(1 + z * z / df) / 2;
  return(out);
}


/*==================================================================================
 * dt_ls - evaluate multivariate t-Student distribution (log scale)
 *
 * args:
 * - x, a point of the support
 * - df, degree of freedom
 * - mu, mean
 * - Sigma, covariance matrix
 *
 * return:
 * - double, log density
 ==================================================================================*/
double dt_ls_mv(arma::vec x,
                double df,
                arma::vec mean,
                arma::mat sigma){
  int d = x.n_elem;
  double out;

  arma::mat rooti  = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double c = lgamma((d + df)/2.0) - lgamma(df/2.0) - (d/2.0) * log(M_PI * df) - 0.5 * arma::det(sigma);

  arma::vec z = rooti * (x - mean);
  out         = c - 0.5 * (df + d)  * std::log1p(arma::sum(z%z) / df);

  return(out);
}
