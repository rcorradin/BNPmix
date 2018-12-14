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
 * Accelerate - UNIVARIATE marginal Polya Urn
 * acceleration step for reshuffling the parameters
 * given an allocation
 *
 * args:
 * - data:  vector of observation
 * - mu:    vector of location component
 * - s2:    vector of scale component
 * - clust: vector of allocation
 * - m0:    mean of location's prior distribution (scalar)
 * - k0:    parameter of NIG on scale of location distribution (scalar)
 * - a0:    shape of prior Gamma distribution on the scale component (scalar)
 * - b0:    rate of prior Gamma distribution on the scale component (scalar)
 *
 * Void function
*/

void accelerate_MAR(arma::vec data,
                    arma::vec &mu,
                    arma::vec &s2,
                    arma::vec clust,
                    double m0,
                    double k0,
                    double a0,
                    double b0){
  for(int j = 0; j < mu.n_elem; j++){
    int nj = sum(clust == j);
    arma::vec tdata = data.elem(arma::find(clust == j));

    double kn = 1.0 / ( (1.0/k0) + nj);
    double mn = kn * ((m0/k0) + sum(tdata));
    double an = a0 + (nj / 2.0);
    double bn = b0 + (((pow(m0, 2)/ k0) + arma::accu(arma::pow(tdata, 2)) - (pow(mn, 2)/ kn)) / 2.0);

    s2[j] = 1.0 / arma::randg(arma::distr_param(an, 1.0 / bn));
    mu[j] = arma::randn() * sqrt(kn * s2[j]) + mn;
  }
}

/*
 * Accelerate - MULTIVARIATE marginal Polya Urn
 * acceleration step for reshuffling the parameters
 * given an allocation
 *
 * args:
 * - data:  matrix of observation
 * - mu:    matrix of location component
 * - s2:    cube of scale component
 * - clust: vector of allocation
 * - m0:    vector of location's prior distribution
 * - k0:    double of NIG on scale of location distribution
 * - S0:    matrix of Inverse Wishart distribution
 * - n0:    degree of freedom of Inverse Wishart distribution
 *
 * Void function
 */

void accelerate_MAR_mv(arma::mat data,
                       arma::mat &mu,
                       arma::cube &s2,
                       arma::vec clust,
                       arma::vec m0,
                       double k0,
                       arma::mat S0,
                       double n0){
  // loop over the different clusters
  for(int j = 0; j < mu.n_rows; j++){
    // initialize itra cluster quantities
    int nj = sum(clust == j);
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
  }
}

/*
 * Clean parameter - UNIVARIATE marginal Polya Urn
 * discard the middle not used values for the clusters and
 * update the correspondent parameters.
 *
 * args:
 * - mu:      vector, each element a mean
 * - s2:      vector, each element a variance
 * - clust:   vector, each (integer) value is the cluster of the corresp. obs.
 *
 * Void function.
*/

void para_clean_MAR(arma::vec &mu,
                    arma::vec &s2,
                    arma::vec &clust) {
  int k = mu.n_elem;

  // for all the used parameters
  for(int i = 0; i < k; i++){

    // if a cluster is empty
    if((int) arma::sum(clust == i) == 0){

      // find the last full cluster, then swap
      for(int j = k; j > i; j--){
        if((int) arma::sum(clust == j) != 0){

          // SWAPPING!!
          clust( arma::find(clust == j) ).fill(i);

          double tmu = mu[i];
          mu[i] = mu[j];
          mu[j] = tmu;

          double ts2 = s2[i];
          s2[i] = s2[j];
          s2[j] = ts2;

          break;
        }
      }
    }
  }

  // reduce dimensions
  int u_bound = 0;
  for(int i = 0; i < k; i++){
    if(arma::accu(clust == i) > 0){
      u_bound += 1;
    }
  }

  // resize object to the correct dimension
  mu.resize(u_bound);
  s2.resize(u_bound);
}

/*
 * Clean parameter - MULTIVARIATE marginal Polya Urn
 * discard the middle not used values for the clusters and
 * update the correspondent parameters.
 *
 * args:
 * - mu:      matrix, each row a mean
 * - s2:      cube, each slice a covariance matrix
 * - clust:   vector, each (integer) value is the cluster of the corresp. obs.
 *
 * Void function.
 */

void para_clean_MAR_mv(arma::mat &mu,
                       arma::cube &s2,
                       arma::vec &clust) {
  int k = mu.n_rows;

  // for all the used parameters
  for(int i = 0; i < k; i++){

    // if a cluster is empty
    if((int) arma::sum(clust == i) == 0){

      // find the last full cluster, then swap
      for(int j = k; j > i; j--){
        if((int) arma::sum(clust == j) != 0){

          // SWAPPING!!
          clust( arma::find(clust == j) ).fill(i);
          mu.swap_rows(i,j);
          s2.slice(i).swap(s2.slice(j));
          break;
        }
      }
    }
  }

  // reduce dimensions
  int u_bound = 0;
  for(int i = 0; i < k; i++){
    if(arma::accu(clust == i) > 0){
      u_bound += 1;
    }
  }

  // resize object to the correct dimension
  mu.resize(u_bound,mu.n_cols);
  s2.resize(s2.n_rows, s2.n_cols, u_bound);
}

/*
 * Update clusters - UNIVARIATE marginal Polya Urn
 *
 * args:
 * - data:        vector of observation
 * - mujoin:      mean values for each component
 * - s2join:      mean values for each component
 * - probjoin:    mean values for each component
 * - clust:       vector of allocation
 *
 * void function
 */

void clust_update_MAR(arma::vec data,
                      arma::vec &mu,
                      arma::vec &s2,
                      arma::vec &clust,
                      double mass,
                      double m0,
                      double k0,
                      double a0,
                      double b0,
                      int iter,
                      arma::vec &new_val){

  // initialize quantities
  int n = clust.n_elem;

  // loop over the observations
  for(int i = 0; i < n; i++){

    bool req_clean = false;
    if(arma::sum(clust == clust[i]) == 1){
      req_clean = true;
    }

    clust(i) = mu.n_elem + 1;
    if(req_clean){
      para_clean_MAR(mu,
                     s2,
                     clust);
    }

    // initialize useful quantities
    int k = mu.n_elem;
    arma::vec probs(k+1);
    probs.fill(0);
    // arma::vec temp_clust = clust;
    // temp_clust(i) = k+1;

    // compute probabilities vector
    for(int j = 0; j < k; j++) {
      int nj = (int) arma::sum(clust == j);
      probs[j] = arma::normpdf(data[i], mu(j), sqrt(s2(j))) * nj;
    }
    probs[k] = mass * dt_ls(data[i], 2 * a0, m0, sqrt(b0 * (1 + k0) / a0));

    // sample new
    int temp_cl = rintnunif(probs);
    clust[i] = temp_cl;

    if(temp_cl == k){
      mu.resize(k+1);
      s2.resize(k+1);

      double kn = 1.0 / ( (1.0/k0) + 1.0);
      double mn = kn * ((m0/k0) + data[i]);
      double an = a0 + (1.0 / 2.0);
      double bn = b0 + (((pow(m0, 2)/ k0) + pow(data[i],2) - (pow(mn, 2)/ kn)) / 2.0);

      s2(k) = 1.0 / arma::randg(arma::distr_param(an, 1.0 / bn));
      mu(k) = arma::randn() * sqrt(kn * s2(k)) + mn;

      new_val(iter)++;
    }
  }
}

/*
 * Update clusters - UNIVARIATE marginal Polya Urn - PY
 *
 * args:
 * - data:        vector of observation
 * - mujoin:      mean values for each component
 * - s2join:      mean values for each component
 * - probjoin:    mean values for each component
 * - clust:       vector of allocation
 * - mass:        DP mass
 * - m0:          mean of location component of base measure
 * - sigma_PY:    second parameter of PY process
 *
 * void function
 */

void clust_update_MAR_PY(arma::vec data,
                         arma::vec &mu,
                         arma::vec &s2,
                         arma::vec &clust,
                         double mass,
                         double m0,
                         double k0,
                         double a0,
                         double b0,
                         int iter,
                         arma::vec &new_val,
                         double sigma_PY){

  // initialize quantities
  int n = clust.n_elem;

  // loop over the observations
  for(int i = 0; i < n; i++){

    bool req_clean = false;
    if(arma::sum(clust == clust[i]) == 1){
      req_clean = true;
    }

    clust(i) = mu.n_elem + 1;
    if(req_clean){
      para_clean_MAR(mu,
                     s2,
                     clust);
    }

    // initialize useful quantities
    int k = mu.n_elem;
    arma::vec probs(k+1);
    probs.fill(0);
    // arma::vec temp_clust = clust;
    // temp_clust(i) = k+1;

    // compute probabilities vector
    for(int j = 0; j < k; j++) {
      int nj = (int) arma::sum(clust == j);
      probs[j] = arma::normpdf(data[i], mu(j), sqrt(s2(j))) * (nj - sigma_PY);
    }
    probs[k] = (mass + k * sigma_PY) * dt_ls(data[i], 2 * a0, m0, sqrt(b0 * (1 + k0) / a0));

    // sample new
    int temp_cl = rintnunif(probs);
    clust[i] = temp_cl;

    if(temp_cl == k){
      mu.resize(k+1);
      s2.resize(k+1);

      double kn = 1.0 / ( (1.0/k0) + 1.0);
      double mn = kn * ((m0/k0) + data[i]);
      double an = a0 + (1.0 / 2.0);
      double bn = b0 + (((pow(m0, 2)/ k0) + pow(data[i],2) - (pow(mn, 2)/ kn)) / 2.0);

      s2(k) = 1.0 / arma::randg(arma::distr_param(an, 1.0 / bn));
      mu(k) = arma::randn() * sqrt(kn * s2(k)) + mn;

      new_val(iter)++;
    }
  }
}

/*
 * Update clusters - MULTIVARIATE marginal Polya Urn
 *
 * args:
 * - data:        vector of observation
 * - mujoin:      mean values for each component
 * - s2join:      mean values for each component
 * - probjoin:    mean values for each component
 * - clust:       vector of allocation
 *
 * void function
 */

void clust_update_MAR_mv(arma::mat data,
                         arma::mat &mu,
                         arma::cube &s2,
                         arma::vec &clust,
                         double mass,
                         arma::vec m0,
                         double k0,
                         arma::mat S0,
                         double n0,
                         int iter,
                         arma::vec &new_val){
  // initialize quantities
  int n = clust.n_elem;
  int d = data.n_cols;

  // loop over the observations
  for(int i = 0; i < n; i++){

    bool req_clean = false;
    if(arma::sum(clust == clust[i]) == 1){
      req_clean = true;
    }

    clust(i) = mu.n_rows + 1;
    if(req_clean){
      para_clean_MAR_mv(mu,
                        s2,
                        clust);
    }

    // initialize useful quantities
    int k = mu.n_rows;
    arma::vec probs(k+1);
    probs.fill(0);

    // compute probabilities vector
    for(int j = 0; j < k; j++) {
      int nj = (int) arma::sum(clust == j);
      arma::mat rooti  = arma::trans(arma::inv(trimatu(arma::chol(s2.slice(j)))));
      arma::vec cdata  = rooti * arma::trans(data.row(i) - mu.row(j)) ;
      double out       = - (d / 2.0) * log(2.0 * M_PI) - 0.5 * arma::sum(cdata%cdata) + arma::sum(log(rooti.diag()));
      probs[j] = exp(out) * nj;
    }
    probs[k] = mass * dt_ls_mv(data.row(i).t(), n0 - d + 1, m0, S0 * (k0 + 1) / (k0 * (n0 - d + 1)));

    // sample new
    int temp_cl = rintnunif(probs);
    clust[i] = temp_cl;

    if(temp_cl == k){
      mu.resize(k+1, mu.n_cols);
      s2.resize(s2.n_rows, s2.n_cols, k+1);

      // update parameters
      double kn = k0 + 1;
      arma::vec mn = ((m0 * k0) + data.row(i).t())/kn;
      double nn = n0 + 1;

      arma::mat Sn = S0 + (data.row(i).t() - mn) * arma::trans(data.row(i).t() - mn) + ((n0+1)/nn) *
        (data.row(i).t() - m0) * arma::trans(data.row(i).t() - m0);

      // sample from the posterior distributions
      s2.slice(k) = arma::inv(arma::wishrnd(arma::inv(Sn), nn));
      mu.row(k)   = arma::trans(arma::mvnrnd(mn, s2.slice(k)/kn));

      new_val[iter]++;
    }
  }

}

/*
 * Update clusters - MULTIVARIATE marginal Polya Urn - PY
 *
 * args:
 * - data:        vector of observation
 * - mujoin:      mean values for each component
 * - s2join:      mean values for each component
 * - probjoin:    mean values for each component
 * - clust:       vector of allocation
 *
 * void function
 */

void clust_update_MAR_PY_mv(arma::mat data,
                            arma::mat &mu,
                            arma::cube &s2,
                            arma::vec &clust,
                            double mass,
                            arma::vec m0,
                            double k0,
                            arma::mat S0,
                            double n0,
                            int iter,
                            arma::vec &new_val,
                            double sigma_PY){
  // initialize quantities
  int n = clust.n_elem;
  int d = data.n_cols;

  // loop over the observations
  for(int i = 0; i < n; i++){

    bool req_clean = false;
    if(arma::sum(clust == clust[i]) == 1){
      req_clean = true;
    }

    clust(i) = mu.n_rows + 1;
    if(req_clean){
      para_clean_MAR_mv(mu,
                        s2,
                        clust);
    }

    // initialize useful quantities
    int k = mu.n_rows;
    arma::vec probs(k+1);
    probs.fill(0);

    // compute probabilities vector
    for(int j = 0; j < k; j++) {
      int nj = (int) arma::sum(clust == j);
      arma::mat rooti  = arma::trans(arma::inv(trimatu(arma::chol(s2.slice(j)))));
      arma::vec cdata  = rooti * arma::trans(data.row(i) - mu.row(j)) ;
      double out       = - (d / 2.0) * log(2.0 * M_PI) - 0.5 * arma::sum(cdata%cdata) + arma::sum(log(rooti.diag()));
      probs[j] = exp(out) * (nj - sigma_PY);
    }
    probs[k] = (mass + k * sigma_PY) * dt_ls_mv(data.row(i).t(), n0 - d + 1, m0, S0 * (k0 + 1) / (k0 * (n0 - d + 1)));

    // sample new
    int temp_cl = rintnunif(probs);
    clust[i] = temp_cl;

    if(temp_cl == k){
      mu.resize(k+1, mu.n_cols);
      s2.resize(s2.n_rows, s2.n_cols, k+1);

      // update parameters
      double kn = k0 + 1;
      arma::vec mn = ((m0 * k0) + data.row(i).t())/kn;
      double nn = n0 + 1;

      arma::mat Sn = S0 + (data.row(i).t() - mn) * arma::trans(data.row(i).t() - mn) + ((n0+1)/nn) *
        (data.row(i).t() - m0) * arma::trans(data.row(i).t() - m0);

      // sample from the posterior distributions
      s2.slice(k) = arma::inv(arma::wishrnd(arma::inv(Sn), nn));
      mu.row(k)   = arma::trans(arma::mvnrnd(mn, s2.slice(k)/kn));

      new_val[iter]++;
    }
  }

}

