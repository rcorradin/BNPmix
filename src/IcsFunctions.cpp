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
#include "CommonUtilities.h"
// [[Rcpp::depends("RcppArmadillo")]]

/*
 * Accelerate - UNIVARIATE conditional Polya Urn Scheme
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

void accelerate_ICS(arma::vec data,
                arma::vec &mu,
                arma::vec &s2,
                arma::vec clust,
                double m0,
                double k0,
                double a0,
                double b0){
  for(arma::uword j = 0; j < mu.n_elem; j++){
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
 * Accelerate - MULTIVARIATE conditional Polya Urn Scheme
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

void accelerate_ICS_mv(arma::mat data,
                       arma::mat &mu,
                       arma::cube &s2,
                       arma::vec clust,
                       arma::vec m0,
                       double k0,
                       arma::mat S0,
                       double n0){
  // loop over the different clusters
  for(arma::uword j = 0; j < mu.n_rows; j++){
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
 * Clean parameter - UNIVARIATE conditional Polya Urn Scheme
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

void para_clean_ICS(arma::vec &mu,
                    arma::vec &s2,
                    arma::vec &clust) {
  int k = mu.n_elem;

  // for all the used parameters
  for(arma::uword i = 0; i < k; i++){

    // if a cluster is empty
    if((int) arma::sum(clust == i) == 0){

      // find the last full cluster, then swap
      for(arma::uword j = k; j > i; j--){
        if((int) arma::sum(clust == j) != 0){

          // swap the correpsonding elements
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
  for(arma::uword i = 0; i < k; i++){
    if(arma::accu(clust == i) > 0){
      u_bound += 1;
    }
  }

  // resize object to the correct dimension
  mu.resize(u_bound);
  s2.resize(u_bound);
}

/*
 * Clean parameter - MULTIVARIATE conditional Polya Urn Scheme
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

void para_clean_ICS_mv(arma::mat &mu,
                       arma::cube &s2,
                       arma::vec &clust) {
  int k = mu.n_rows;

  // for all the used parameters
  for(arma::uword i = 0; i < k; i++){

    // if a cluster is empty
    if((int) arma::sum(clust == i) == 0){

      // find the last full cluster, then swap
      for(arma::uword j = k; j > i; j--){
        if((int) arma::sum(clust == j) != 0){

          // swap the corresponding elements
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
  for(arma::uword i = 0; i < k; i++){
    if(arma::accu(clust == i) > 0){
      u_bound += 1;
    }
  }

  // resize object to the correct dimension
  mu.resize(u_bound,mu.n_cols);
  s2.resize(s2.n_rows, s2.n_cols, u_bound);
}

/*
 * Simulate finite distribution - UNIVARIATE - sequence from Pitman Yor process
 *
 * args:
 * - mutemp:      mean values for each temp component
 * - s2temp:      variance values for temp each component
 * - freqtemp:    frequency values for temp each component
 * - mass:        mass of Dirichlet process
 * - m0:          mean of location's prior distribution (scalar)
 * - k0:          parameter of NIG on scale of location distribution (scalar)
 * - a0:          shape of prior Gamma distribution on the scale component (scalar)
 * - b0:          rate of prior Gamma distribution on the scale component (scalar)
 * - napprox:     number of approximating values
 * - sigma_PY:    second parameter of Pitman-Yor process
 *
 * Void function.
 */

void simu_trunc_PY(arma::vec &mutemp,
                    arma::vec &s2temp,
                    arma::vec &freqtemp,
                    double mass,
                    double m0,
                    double k0,
                    double a0,
                    double b0,
                    int napprox,
                    double sigma_PY){

  // By independence of atoms and jumps, sample
  // the jumps, then the entire vector of atoms
  // resize the objects to size 1
  freqtemp.resize(1);

  // initialize the first element
  freqtemp.fill(1);
  int k = 1;

  // generate napprox values with ties
  for(arma::uword j = 1; j < napprox; j++){

    int temp_cl = rintnunifw(freqtemp - sigma_PY, mass + freqtemp.n_elem * sigma_PY);

    if(temp_cl < (k - 1)){

      // if is an existing one, increase the freq
      freqtemp[temp_cl] += 1;

    } else {

      // if is a new one, generate the new parameters
      freqtemp.resize(k + 1);
      freqtemp[k] = 1;
      k += 1;

    }
  }

  mutemp.resize(k);
  s2temp.resize(k);
  s2temp = 1.0 / arma::randg(k, arma::distr_param(a0, 1.0 / b0));
  mutemp = arma::randn(k) % sqrt(k0 * s2temp) + m0;

}

/*
 * Simulate finite distribution - MULTIVARIATE - sequence from Pitman Yor process
 *
 * args:
 * - mutemp:      mean values for each temp component
 * - s2temp:      variance values for temp each component
 * - freqtemp:    frequency values for temp each component
 * - mass:        mass of Dirichlet process
 * - m0:          vector of location's prior distribution
 * - k0:          double of NIG on scale of location distribution
 * - S0:          matrix of Inverse Wishart distribution
 * - n0:          degree of freedom of Inverse Wishart distribution
 * - napprox:     number of approximating values
 * - sigma_PY:    second parameter of Pitman-Yor process
 *
 * Void function.
 */

void simu_trunc_PY_mv(arma::mat &mutemp,
                       arma::cube &s2temp,
                       arma::vec &freqtemp,
                       double mass,
                       arma::vec m0,
                       double k0,
                       arma::mat S0,
                       double n0,
                       int napprox,
                       double sigma_PY){

  // resize the objects to size 1
  freqtemp.resize(1);

  // initialize the first element
  freqtemp.fill(1);
  int k = 1;

  // generate napprox values with ties
  for(arma::uword j = 1; j < napprox; j++){
    int temp_cl = rintnunifw(freqtemp - sigma_PY, mass + freqtemp.n_elem * sigma_PY);
    if(temp_cl < (k - 1)){

      // if is an existing one, increase the freq
      freqtemp[temp_cl] += 1;

    } else {

      // if is a new one, generate the new parameters
      freqtemp.resize(k + 1);
      freqtemp[k] = 1;
      k += 1;

    }
  }

  mutemp.resize(k, mutemp.n_cols);
  s2temp.resize(s2temp.n_rows, s2temp.n_cols, k);
  for(arma::uword j = 0; j < k; j++){
    s2temp.slice(j) = arma::inv(arma::wishrnd(arma::inv(S0), n0));
    mutemp.row(j) = arma::trans(arma::mvnrnd(m0, s2temp.slice(j)/k0));
  }
}


/*
 * Update clusters - UNIVARIATE conditional Polya Urn Scheme
 *
 * args:
 * - data:        vector of observation
 * - mujoin:      mean values for each component
 * - s2join:      mean values for each component
 * - probjoin:    mean values for each component
 * - clust:       vector of allocation
 * - max_val:     vector, number of already existent atoms
 * - iter:        current iteration
 * - new_val:     vector of new values
 *
 * void function
 */

void clust_update_ICS(arma::vec data,
                      arma::vec mujoin,
                      arma::vec s2join,
                      arma::vec probjoin,
                      arma::vec &clust,
                      int max_val,
                      int iter,
                      arma::vec &new_val){

  // initialize the quantities
  int n = clust.n_elem;
  int k = probjoin.n_elem;
  arma::vec probs_upd(k);

  // loop over the observations
  for(arma::uword i = 0; i < n; i++){

    // loop over the components
    for(arma::uword j = 0; j < k; j++){
      probs_upd[j] = probjoin[j] * arma::normpdf(data[i], mujoin[j], sqrt(s2join[j]));
    }

    // sample the allocation for the current observation
    clust[i] = rintnunif(probs_upd);
    if(clust[i] > max_val - 1){
      new_val[iter]++;
    }
  }
}

/*
 * Update clusters - MULTIVARIATE conditional Polya Urn Scheme
 *
 * args:
 * - data:        vector of observation
 * - mujoin:      mean values for each component
 * - s2join:      mean values for each component
 * - probjoin:    mean values for each component
 * - clust:       vector of allocation
 * - max_val:     vector, number of already existent atoms
 * - iter:        current iteration
 * - new_val:     vector of new values
 *
 * void function
 */

void clust_update_ICS_mv(arma::mat data,
                         arma::mat mujoin,
                         arma::cube s2join,
                         arma::vec probjoin,
                         arma::vec &clust,
                         int max_val,
                         int iter,
                         arma::vec &new_val){

  // initialize the quantities
  double d = (double) data.n_cols;
  int n = clust.n_elem;
  int k = probjoin.n_elem;
  arma::vec probs_upd(k);

  // loop over the observations
  for(arma::uword i = 0; i < n; i++){

    // loop over the components
    for(arma::uword j = 0; j < k; j++){
      arma::mat rooti  = arma::trans(arma::inv(trimatu(arma::chol(s2join.slice(j)))));
      arma::vec cdata  = rooti * arma::trans(data.row(i) - mujoin.row(j)) ;
      double out       = - (d / 2.0) * log(2.0 * M_PI) - 0.5 * arma::sum(cdata%cdata) +
                            arma::sum(log(rooti.diag()));
      probs_upd[j]     = probjoin[j] * exp(out);
    }

    // sample the allocation for the current observations
    clust[i] = rintnunif(probs_upd);
    if(clust[i] > max_val - 1){
      new_val[iter]++;
    }
  }
}

