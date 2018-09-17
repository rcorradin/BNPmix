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

// NON-Uniform discrete distribution
int rintnunif(arma::vec weights){
  double u = arma::randu();
  arma::vec probs = weights / sum(weights);
  probs = arma::cumsum(probs);
  
  for(int k = 0; k < probs.n_elem; k++) {
    if(u <= probs[k]) {
      return k;
    }
  }
}

// normalized freq + mass 
int rintnunifw(arma::vec freq,
               double mass){
  arma::vec weights(freq);
  weights.resize(freq.n_elem + 1);
  weights[freq.n_elem] = mass;
  arma::vec probs = arma::cumsum(weights) / sum(weights);
  double u = arma::randu();
  
  for(int k = 0; k < probs.n_elem; k++) {
    if(u <= probs[k]) {
      return k;
    }
  }
}

// Dirichlet distribution
arma::vec rdirich_mass(arma::vec freq, 
                       double mass){
  
  arma::vec weights(freq);
  weights.resize(freq.n_elem + 1);
  weights[freq.n_elem] = mass;
  arma::vec result(weights.n_elem);

  for(int j = 0; j < weights.n_elem; j++){
    result[j] = arma::randg(1, arma::distr_param(weights[j], 1.0))[0];
  }
  return(result / sum(result));
}

// univariate t student density
double dt_ls(double x,
             double df,
             double mu,
             double sigma){
  double z = (x - mu)/sigma;
  double out = lgamma((df + 1) / 2) - log(sqrt(M_PI * df)) - 
    lgamma(df / 2) - (df + 1) * log(1 + z * z / df) / 2;
  return(exp(out));
}

// multivariate t distribution density 
double dt_ls_mv(arma::vec x,
                double df,
                arma::vec mean,
                arma::mat sigma){
  int d = x.n_elem;
  double out;
  
  arma::mat rooti  = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum  = arma::sum(log(rooti.diag()));
  double c = lgamma((d + df)/2.0) - lgamma(df/2.0) - (d/2.0) * log(M_PI * df) + rootisum;
  
  arma::vec z = rooti * (x - mean) ;
  out         = c - 0.5 * (df + d)  * std::log1p(arma::sum(z%z) / df);
  
  return(exp(out));
}