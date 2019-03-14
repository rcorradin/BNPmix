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
#include "Distributions.h"
// [[Rcpp::depends("RcppArmadillo")]]

/*
 * Accelerate - UNIVARIATE conditional Slice sampler
 * acceleration step for reshuffling the parameters
 * given an allocation
 *
 * args:
 * - data:  vector of observation
 * - mu:    vector of location component
 * - s2:    vector of scale component
 * - v:     vector of stick break components
 * - w:     vector of stick weights
 * - clust: vector of allocation
 * - m0:    mean of location's prior distribution (scalar)
 * - k0:    parameter of NIG on scale of location distribution (scalar)
 * - a0:    shape of prior Gamma distribution on the scale component (scalar)
 * - b0:    rate of prior Gamma distribution on the scale component (scalar)
 * - mass:  mass of Dirichlet process
 *
 * Void function
 */

void accelerate_SLI(arma::vec data,
                   arma::vec &mu,
                   arma::vec &s2,
                   arma::vec &v,
                   arma::vec &w,
                   arma::vec clust,
                   double m0,
                   double k0,
                   double a0,
                   double b0,
                   double mass){
  double xtemp;
  double ytemp;

  for(arma::uword j = 0; j < mu.n_elem; j++){
    int nj  = arma::accu(clust == j);
    int ngj = arma::accu(clust > j);
    arma::vec tdata = data.elem(arma::find(clust == j));

    double kn = 1.0 / ( (1.0/k0) + nj);
    double mn = kn * ((m0/k0) + sum(tdata));
    double an = a0 + (nj / 2.0);
    double bn = b0 + (((pow(m0, 2)/ k0) + arma::accu(arma::pow(tdata, 2)) - (pow(mn, 2)/ kn)) / 2.0);

    s2[j] = 1.0 / arma::randg(arma::distr_param(an, 1.0 / bn));
    mu[j] = arma::randn() * sqrt(kn * s2[j]) + mn;
    xtemp = arma::randg(arma::distr_param(1.0 + nj, 1.0));
    ytemp = arma::randg(arma::distr_param(mass + ngj, 1.0));
    v[j]  = xtemp / (xtemp + ytemp);

    if(j == 0){
      w[j] = v[j];
    }else{
      w[j] = v[j] * ((1 - v[j-1]) * w[j - 1]) / v[j-1];
    }
  }
}

/*
 * Accelerate - UNIVARIATE conditional Slice sampler - PY
 * acceleration step for reshuffling the parameters
 * given an allocation
 *
 * args:
 * - data:  vector of observation
 * - mu:    vector of location component
 * - s2:    vector of scale component
 * - v:     vector of stick break components
 * - w:     vector of stick weights
 * - clust: vector of allocation
 * - m0:    mean of location's prior distribution (scalar)
 * - k0:    parameter of NIG on scale of location distribution (scalar)
 * - a0:    shape of prior Gamma distribution on the scale component (scalar)
 * - b0:    rate of prior Gamma distribution on the scale component (scalar)
 * - mass:  mass of PY process
 * - sigma_PY: second parameter PY
 *
 * Void function
 */

void accelerate_SLI_PY(arma::vec data,
                       arma::vec &mu,
                       arma::vec &s2,
                       arma::vec &v,
                       arma::vec &w,
                       arma::vec clust,
                       double m0,
                       double k0,
                       double a0,
                       double b0,
                       double mass,
                       double sigma_PY){
  double xtemp;
  double ytemp;

  for(arma::uword j = 0; j < mu.n_elem; j++){
    int nj  = arma::accu(clust == j);
    int ngj = arma::accu(clust > j);
    arma::vec tdata = data.elem(arma::find(clust == j));

    double kn = 1.0 / ( (1.0/k0) + nj);
    double mn = kn * ((m0/k0) + sum(tdata));
    double an = a0 + (nj / 2.0);
    double bn = b0 + (((pow(m0, 2)/ k0) + arma::accu(arma::pow(tdata, 2)) - (pow(mn, 2)/ kn)) / 2.0);

    s2[j] = 1.0 / arma::randg(arma::distr_param(an, 1.0 / bn));
    mu[j] = arma::randn() * sqrt(kn * s2[j]) + mn;
    xtemp = arma::randg(arma::distr_param(1.0 - sigma_PY + nj, 1.0));
    ytemp = arma::randg(arma::distr_param(mass + (j + 1) * sigma_PY + ngj, 1.0));
    v[j]  = xtemp / (xtemp + ytemp);

    if(j == 0){
      w[j] = v[j];
    }else{
      w[j] = v[j] * ((1 - v[j-1]) * w[j - 1]) / v[j-1];
    }
  }
}

/*
 * Accelerate - MULTIVARIATE conditional Slice sampler
 * acceleration step for reshuffling the parameters
 * given an allocation
 *
 * args:
 * - data:  matrix of observation
 * - mu:    matrix of location component
 * - s2:    cube of scale component
 * - v:     vector of stick break components
 * - w:     vector of stick weights
 * - clust: vector of allocation
 * - m0:    vector of location's prior distribution
 * - k0:    double of NIG on scale of location distribution
 * - S0:    matrix of Inverse Wishart distribution
 * - n0:    degree of freedom of Inverse Wishart distribution
 * - mass:  mass of Dirichlet process
 *
 * Void function
 */

void accelerate_SLI_mv(arma::mat data,
                       arma::mat &mu,
                       arma::cube &s2,
                       arma::vec &v,
                       arma::vec &w,
                       arma::vec clust,
                       arma::vec m0,
                       double k0,
                       arma::mat S0,
                       double n0,
                       double mass){
  double xtemp;
  double ytemp;

  // loop over the clusters
  for(arma::uword j = 0; j < mu.n_rows; j++){
    // initialize intra cluster quantities
    int nj  = arma::accu(clust == j);
    int ngj = arma::accu(clust > j);
    arma::mat tdata = data.rows(arma::find(clust == j));

    // update parameters
    double kn = k0 + nj;
    arma::vec mn = ((m0 * k0) + arma::trans(sum(tdata, 0)))/kn;
    arma::mat cdata = tdata - arma::repmat(mean(tdata, 0), nj, 1);
    double nn = n0 + nj;
    arma::mat Sn = S0 + arma::trans(cdata) * cdata + ((n0+nj)/nn) *
      (arma::trans(mean(tdata, 0)) - m0) * arma::trans(arma::trans(mean(tdata, 0)) - m0);

    // sample from the posterior distributions
    s2.slice(j) = arma::inv(arma::wishrnd(arma::inv(Sn), nn));
    mu.row(j)   = arma::trans(arma::mvnrnd(mn, s2.slice(j)/kn));
    xtemp = arma::randg(arma::distr_param(1.0 + nj, 1.0));
    ytemp = arma::randg(arma::distr_param(mass + ngj, 1.0));
    v[j]  = xtemp / (xtemp + ytemp);

    if(j == 0){
      w[j] = v[j];
    }else{
      w[j] = v[j] * ((1 - v[j-1]) * w[j - 1]) / v[j-1];
    }
  }
}

/*
 * Accelerate - MULTIVARIATE conditional Slice sampler - PY
 * acceleration step for reshuffling the parameters
 * given an allocation
 *
 * args:
 * - data:  matrix of observation
 * - mu:    matrix of location component
 * - s2:    cube of scale component
 * - v:     vector of stick break components
 * - w:     vector of stick weights
 * - clust: vector of allocation
 * - m0:    vector of location's prior distribution
 * - k0:    double of NIG on scale of location distribution
 * - S0:    matrix of Inverse Wishart distribution
 * - n0:    degree of freedom of Inverse Wishart distribution
 * - mass:  mass of Dirichlet process
 * - sigma_PY: second parameter PY
 *
 * Void function
 */

void accelerate_SLI_PY_mv(arma::mat data,
                          arma::mat &mu,
                          arma::cube &s2,
                          arma::vec &v,
                          arma::vec &w,
                          arma::vec clust,
                          arma::vec m0,
                          double k0,
                          arma::mat S0,
                          double n0,
                          double mass,
                          double sigma_PY){
  double xtemp;
  double ytemp;

  // loop over the clusters
  for(arma::uword j = 0; j < mu.n_rows; j++){
    // initialize intra cluster quantities
    int nj  = arma::accu(clust == j);
    int ngj = arma::accu(clust > j);
    arma::mat tdata = data.rows(arma::find(clust == j));

    // update parameters
    double kn = k0 + nj;
    arma::vec mn = ((m0 * k0) + arma::trans(sum(tdata, 0)))/kn;
    arma::mat cdata = tdata - arma::repmat(mean(tdata, 0), nj, 1);
    double nn = n0 + nj;
    arma::mat Sn = S0 + arma::trans(cdata) * cdata + ((n0+nj)/nn) *
      (arma::trans(mean(tdata, 0)) - m0) * arma::trans(arma::trans(mean(tdata, 0)) - m0);

    // sample from the posterior distributions
    s2.slice(j) = arma::inv(arma::wishrnd(arma::inv(Sn), nn));
    mu.row(j)   = arma::trans(arma::mvnrnd(mn, s2.slice(j)/kn));
    xtemp = arma::randg(arma::distr_param(1.0 - sigma_PY + nj, 1.0));
    ytemp = arma::randg(arma::distr_param(mass + (j + 1) * sigma_PY + ngj, 1.0));
    v[j]  = xtemp / (xtemp + ytemp);

    if(j == 0){
      w[j] = v[j];
    }else{
      w[j] = v[j] * ((1 - v[j-1]) * w[j - 1]) / v[j-1];
    }
  }
}

/*
 * Clean parameter - UNIVARIATE conditional Slice sampler
 * discard the middle not used values for the clusters and
 * update the correspondent parameters.
 *
 * args:
 * - mu:      vector, each element a mean
 * - s2:      vector, each element a variance
 * - clust:   vector, each (integer) value is the cluster of the corresp. obs.
 * - v:       vector of stick break components
 * - w:       vector of stick weights
 *
 * Void function.
 */

void para_clean_SLI(arma::vec &mu,
                    arma::vec &s2,
                    arma::vec &clust,
                    arma::vec &v,
                    arma::vec &w) {
  int k = mu.n_elem;

  // for all the used parameters
  for(arma::uword i = 0; i < k; i++){

    // if a cluster is empty
    if((int) arma::sum(clust == i) == 0){

      // find the last full cluster, then swap
      for(arma::uword j = k; j > i; j--){
        if((int) arma::sum(clust == j) != 0){

          // swap the corresponding elements
          clust( arma::find(clust == j) ).fill(i);

          double tmu = mu[i];
          mu[i] = mu[j];
          mu[j] = tmu;

          double ts2 = s2[i];
          s2[i] = s2[j];
          s2[j] = ts2;

          double tv = v[i];
          v[i] = v[j];
          v[j] = tv;

          double tw = w[i];
          w[i] = w[j];
          w[j] = tw;

          break;
        }
      }
    }
  }

  // reduce dimensions
  int u_bound = 0;
  for(arma::uword i = 0; i < k; i++){
    if(arma::accu(clust == i) > 0){
      u_bound += 1;
    }
  }

  mu.resize(u_bound);
  s2.resize(u_bound);
  v.resize(u_bound);
  w.resize(u_bound);
}


/*
 * Clean parameter - MULTIVARIATE conditional Slice sampler
 * discard the middle not used values for the clusters and
 * update the correspondent parameters.
 *
 * args:
 * - mu:      matrix, each row a mean
 * - s2:      cube, each slice a covariance matrix
 * - clust:   vector, each (integer) value is the cluster of the corresp. obs.
 * - v:       vector of stick break components
 * - w:       vector of stick weights
 *
 * Void function.
 */

void para_clean_SLI_mv(arma::mat &mu,
                       arma::cube &s2,
                       arma::vec &clust,
                       arma::vec &v,
                       arma::vec &w) {
  int k = mu.n_rows;

  // for all the used parameters
  for(arma::uword i = 0; i < k; i++){

    // if a cluster is empty
    if((int) arma::sum(clust == i) == 0){

      // find the last full cluster, then swap
      for(arma::uword j = k; j > i; j--){
        if((int) arma::sum(clust == j) != 0){

          clust( arma::find(clust == j) ).fill(i);

          // swap the corresponding elements
          clust( arma::find(clust == j) ).fill(i);
          mu.swap_rows(i,j);
          s2.slice(i).swap(s2.slice(j));

          double tv = v[i];
          v[i] = v[j];
          v[j] = tv;

          double tw = w[i];
          w[i] = w[j];
          w[j] = tw;

          break;
        }
      }
    }
  }

  // reduce dimensions
  int u_bound = 0;
  for(arma::uword i = 0; i < k; i++){
    if(arma::accu(clust == i) > 0){
      u_bound += 1;
    }
  }

  mu.resize(u_bound, mu.n_cols);
  s2.resize(s2.n_rows, s2.n_cols, u_bound);
  v.resize(u_bound);
  w.resize(u_bound);
}


/*
 * grow parameters - UNIVARIATE conditional Slice sampler
 * growing up the parameter vectors
 * till reaching the condition sum(w) > u_i, for all i
 *
 * args:
 * - mu:      vector, each element a mean
 * - s2:      vector, each element a variance
 * - v:       vector of stick break components
 * - w:       vector of stick weights
 * - u:       vector of uniform values
 * - m0:      mean of location's prior distribution (scalar)
 * - k0:      parameter of NIG on scale of location distribution (scalar)
 * - a0:      shape of prior Gamma distribution on the scale component (scalar)
 * - b0:      rate of prior Gamma distribution on the scale component (scalar)
 * - mass:    mass of Dirichlet process
 * - n:       number of observations
 *
 * Void function
 */

void grow_param_SLI(arma::vec &mu,
                    arma::vec &s2,
                    arma::vec &v,
                    arma::vec &w,
                    arma::vec u,
                    double m0,
                    double k0,
                    double a0,
                    double b0,
                    double mass,
                    int n){
  double xtemp;
  double ytemp;
  double w_sum = arma::accu(w);
  int k_old = mu.n_elem;

  while(sum(1 - u < w_sum) < n){

    int k = w.n_elem;
    v.resize(k+1);
    w.resize(k+1);

    xtemp = arma::randg(arma::distr_param(1.0, 1.0));
    ytemp = arma::randg(arma::distr_param(mass, 1.0));
    v[k]  = xtemp / (xtemp + ytemp);

    if(k == 0){
      w[k] = v[k];
    }else{
      w[k] = v[k] * ((1 - v[k-1]) * w[k - 1]) / v[k-1];
    }
    w_sum = arma::accu(w);
  }

  if(w.n_elem > k_old){
    int new_val = w.n_elem - k_old;

    arma::vec s2_temp = 1.0 / arma::randg(new_val, arma::distr_param(a0, 1.0 / b0));
    arma::vec mu_temp = arma::randn(new_val) % sqrt(k0 * s2_temp) + m0;

    mu = arma::join_cols(mu, mu_temp);
    s2 = arma::join_cols(s2, s2_temp);
  }
}

/*
 * grow parameters - UNIVARIATE conditional Slice sampler - PY
 * growing up the parameter vectors
 * till reaching the condition sum(w) > u_i, for all i
 *
 * args:
 * - mu:      vector, each element a mean
 * - s2:      vector, each element a variance
 * - v:       vector of stick break components
 * - w:       vector of stick weights
 * - u:       vector of uniform values
 * - m0:      mean of location's prior distribution (scalar)
 * - k0:      parameter of NIG on scale of location distribution (scalar)
 * - a0:      shape of prior Gamma distribution on the scale component (scalar)
 * - b0:      rate of prior Gamma distribution on the scale component (scalar)
 * - mass:    mass of Dirichlet process
 * - n:       number of observations
 * - sigma_PY: second parameter PY
 *
 * Void function
 */

void grow_param_SLI_PY(arma::vec &mu,
                       arma::vec &s2,
                       arma::vec &v,
                       arma::vec &w,
                       arma::vec u,
                       double m0,
                       double k0,
                       double a0,
                       double b0,
                       double mass,
                       int n,
                       double sigma_PY){
  double xtemp;
  double ytemp;
  double w_sum = arma::accu(w);
  int k_old = mu.n_elem;

  while(sum(1 - u < w_sum) < n){

    int k = w.n_elem;
    v.resize(k + 1);
    w.resize(k + 1);

    xtemp = arma::randg(arma::distr_param(1.0 - sigma_PY, 1.0));
    ytemp = arma::randg(arma::distr_param(mass + (k + 1) * sigma_PY, 1.0));
    v[k]  = xtemp / (xtemp + ytemp);

    if(k == 0){
      w[k] = v[k];
    }else{
      w[k] = v[k] * ((1 - v[k-1]) * w[k - 1]) / v[k-1];
    }
    w_sum = arma::accu(w);
  }

  if(w.n_elem > k_old){
    int new_val = w.n_elem - k_old;

    arma::vec s2_temp = 1.0 / arma::randg(new_val, arma::distr_param(a0, 1.0 / b0));
    arma::vec mu_temp = arma::randn(new_val) % sqrt(k0 * s2_temp) + m0;

    mu = arma::join_cols(mu, mu_temp);
    s2 = arma::join_cols(s2, s2_temp);
  }
}

/*
 * grow parameters - MULTIVARIATE conditional Slice sampler
 * growing up the parameter vectors
 * till reaching the condition sum(w) > u_i, for all i
 *
 * args:
 * - mu:      matrix, each row a mean
 * - s2:      cube, each slice a covariance matrix
 * - v:       vector of stick break components
 * - w:       vector of stick weights
 * - u:       vector of uniform values
 * - m0:      vector of location's prior distribution
 * - k0:      double of NIG on scale of location distribution
 * - S0:      matrix of Inverse Wishart distribution
 * - n0:      degree of freedom of Inverse Wishart distribution
 * - mass:    mass of Dirichlet process
 * - n:       number of observations
 *
 * Void function
 */

void grow_param_SLI_mv(arma::mat &mu,
                       arma::cube &s2,
                       arma::vec &v,
                       arma::vec &w,
                       arma::vec u,
                       arma::vec m0,
                       double k0,
                       arma::mat S0,
                       double n0,
                       double mass,
                       int n){
  double xtemp;
  double ytemp;
  double w_sum = arma::accu(w);
  int k_old = mu.n_rows;

  while(sum(1 - u < w_sum) < n){

    int k = w.n_elem;
    v.resize(k+1);
    w.resize(k+1);

    xtemp = arma::randg(arma::distr_param(1.0, 1.0));
    ytemp = arma::randg(arma::distr_param(mass, 1.0));
    v[k]  = xtemp / (xtemp + ytemp);

    if(k == 0){
      w[k] = v[k];
    }else{
      w[k] = v[k] * ((1 - v[k-1]) * w[k - 1]) / v[k-1];
    }
    w_sum = arma::accu(w);
  }

  int k_new = w.n_elem;

  mu.resize(k_new, mu.n_cols);
  s2.resize(s2.n_rows, s2.n_cols, k_new);

  for(arma::uword j = k_old; j < k_new; j++){
    s2.slice(j) = arma::inv(arma::wishrnd(arma::inv(S0), n0));
    mu.row(j) = arma::trans(arma::mvnrnd(m0, s2.slice(j)/k0));
  }
}

/*
 * grow parameters - MULTIVARIATE conditional Slice sampler - PY
 * growing up the parameter vectors
 * till reaching the condition sum(w) > u_i, for all i
 *
 * args:
 * - mu:      matrix, each row a mean
 * - s2:      cube, each slice a covariance matrix
 * - v:       vector of stick break components
 * - w:       vector of stick weights
 * - u:       vector of uniform values
 * - m0:      vector of location's prior distribution
 * - k0:      double of NIG on scale of location distribution
 * - S0:      matrix of Inverse Wishart distribution
 * - n0:      degree of freedom of Inverse Wishart distribution
 * - mass:    mass of Dirichlet process
 * - n:       number of observations
 * - sigma_PY: second parameter PY
 *
 * Void function
 */

void grow_param_SLI_PY_mv(arma::mat &mu,
                          arma::cube &s2,
                          arma::vec &v,
                          arma::vec &w,
                          arma::vec u,
                          arma::vec m0,
                          double k0,
                          arma::mat S0,
                          double n0,
                          double mass,
                          int n,
                          double sigma_PY){
  double xtemp;
  double ytemp;
  double w_sum = arma::accu(w);
  int k_old = mu.n_rows;

  while(sum(1 - u < w_sum) < n){

    int k = w.n_elem;
    v.resize(k+1);
    w.resize(k+1);

    xtemp = arma::randg(arma::distr_param(1.0 - sigma_PY, 1.0));
    ytemp = arma::randg(arma::distr_param(mass + (k + 1) * sigma_PY, 1.0));
    v[k]  = xtemp / (xtemp + ytemp);

    if(k == 0){
      w[k] = v[k];
    }else{
      w[k] = v[k] * ((1 - v[k-1]) * w[k - 1]) / v[k-1];
    }
    w_sum = arma::accu(w);
  }

  int k_new = w.n_elem;

  mu.resize(k_new, mu.n_cols);
  s2.resize(s2.n_rows, s2.n_cols, k_new);

  for(arma::uword j = k_old; j < k_new; j++){
    s2.slice(j) = arma::inv(arma::wishrnd(arma::inv(S0), n0));
    mu.row(j) = arma::trans(arma::mvnrnd(m0, s2.slice(j)/k0));
  }
}

/*
 * update u - conditional Slice sampler
 * growing up the parameter vectors
 * till reaching the condition sum(w) > u_i, for all i
 *
 * args:
 * - clust:   vector, each (integer) value is the cluster of the corresp. obs.
 * - w:       vector of stick weights
 * - u:       vector of uniform values
 *
 * Void function
 */

void update_u_SLI(arma::vec clust,
                  arma::vec w,
                  arma::vec &u){
  int nel = clust.n_elem;

  for(arma::uword el = 0; el < nel; el++){
    u[el] = arma::randu() * w(clust[el]);
  }
}

/*
 * update cluster - UNIVARIATE conditional Slice sampler
 *
 * args:
 * - data:  vector of observation
 * - mu:      vector, each element a mean
 * - s2:      vector, each element a variance
 * - clust:   vector, each (integer) value is the cluster of the corresp. obs.
 * - w:       vector of stick weights
 * - u:       vector of uniform values
 *
 * Void function
 */

void update_cluster_SLI(arma::vec data,
                        arma::vec mu,
                        arma::vec s2,
                        arma::vec &clust,
                        arma::vec w,
                        arma::vec u,
                        int max_val,
                        int iter,
                        arma::vec &new_val){
  int n = data.n_elem;
  int k = mu.n_elem;
  arma::uvec index_use;
  arma::uvec index = arma::regspace<arma::uvec>(0, k - 1);
  arma::vec probs;
  int siz;
  int sampled;

  for(arma::uword i = 0; i < n; i++){
    siz = 0;
    index_use.resize(1);
    for(arma::uword r = 0; r < k; r++){
      if(w[r] > u[i]){
        siz++;
        index_use.resize(siz);
        index_use[siz - 1] = index[r];
      }
    }

    if(index_use.n_elem == 1){
      clust[i] = index_use[0];
      if(clust[i] > max_val - 1){
        new_val[iter]++;
      }
    } else {
      probs.resize(index_use.n_elem);
      for(arma::uword j = 0; j < index_use.n_elem; j++){
        probs[j] = arma::normpdf(data[i], mu(index_use[j]), sqrt(s2(index_use[j])));
      }
      sampled = rintnunif(probs);
      clust[i] = index_use[sampled];
      if(clust[i] > max_val - 1){
        new_val[iter]++;
      }
    }
  }
}

/*
 * update cluster - MULTIVARIATE conditional Slice sampler
 *
 * args:
 * - data:  matrix of observation
 * - mu:    matrix of location component
 * - s2:    cube of scale component
 * - clust:   vector, each (integer) value is the cluster of the corresp. obs.
 * - w:       vector of stick weights
 * - u:       vector of uniform values
 *
 * Void function
 */

void update_cluster_SLI_mv(arma::mat data,
                           arma::mat mu,
                           arma::cube s2,
                           arma::vec &clust,
                           arma::vec w,
                           arma::vec u,
                           int max_val,
                           int iter,
                           arma::vec &new_val){
  int n = data.n_rows;
  int k = mu.n_rows;
  int d = data.n_cols;
  arma::uvec index_use;
  arma::uvec index = arma::regspace<arma::uvec>(0, k - 1);
  arma::vec probs;
  int siz;
  int sampled;

  for(arma::uword i = 0; i < n; i++){
    siz = 0;
    index_use.resize(1);
    for(arma::uword r = 0; r < k; r++){
      if(w[r] > u[i]){
        siz++;
        index_use.resize(siz);
        index_use[siz - 1] = index[r];
      }
    }

    if(index_use.n_elem == 1){
      clust[i] = index_use[0];
      if(clust[i] > max_val - 1){
        new_val[iter]++;
      }
    } else {
      probs.resize(index_use.n_elem);
      for(arma::uword j = 0; j < index_use.n_elem; j++){
        arma::mat rooti  = arma::trans(arma::inv(trimatu(arma::chol(s2.slice(index_use[j])))));
        arma::vec cdata  = rooti * arma::trans(data.row(i) - mu.row(index_use[j])) ;
        probs[j]         = exp(- (d / 2.0) * log(2.0 * M_PI) - 0.5 * arma::sum(cdata%cdata) + arma::sum(log(rooti.diag())));
      }
      sampled = rintnunif(probs);
      clust[i] = index_use[sampled];
      if(clust[i] > max_val - 1){
        new_val[iter]++;
      }
    }
  }
}
