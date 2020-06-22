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
#include "Distributions.h"
// [[Rcpp::depends("RcppArmadillo")]]

/*----------------------------------------------------------------------
 *
 * UNIVARIATE
 * LOCATION KERNEL
 * marginal functions
 *
 *----------------------------------------------------------------------
 */

/*==================================================================================
 * Accelerate - UNIVARIATE marginal sampler - LOCATION
 * acceleration step for reshuffling the parameters
 * given an allocation
 *
 * args:
 * - data:  vector of observation
 * - mu:    vector of location component
 * - s2:    scale component
 * - clust: vector of allocation
 * - m0:    mean of location's prior distribution (scalar)
 * - s20:   prior variance of the location component
 * - a0:    shape of prior Gamma distribution on the scale component (scalar)
 * - b0:    rate of prior Gamma distribution on the scale component (scalar)
 *
 * Void function
 ==================================================================================*/

void accelerate_MAR_L(arma::vec data,
                      arma::vec &mu,
                      double s2,
                      arma::vec clust,
                      double m0,
                      double s20,
                      double a0,
                      double b0){
  double s_temp, m_temp;
  double accu_temp = 0.0;
  int nj;

  for(arma::uword j = 0; j < mu.n_elem; j++){
    nj = sum(clust == j);
    arma::vec tdata = data.elem(arma::find(clust == j));

    s_temp = 1 / (nj/s2 + 1/s20);
    m_temp = s_temp * (m0 / s20 + arma::accu(tdata)/s2);

    mu[j] = arma::randn() * sqrt(s_temp) + m_temp;
    accu_temp += arma::accu(pow(tdata - mu[j], 2));
  }

  s2 = 1.0 / arma::randg(arma::distr_param((a0 + data.n_elem / 2), 1.0 / (b0 + accu_temp / 2)));
}

/*==================================================================================
 * Hyper-accelerate - UNIVARIATE marginal sampler - LOCATION
 * hyper-acceleration step for reshuffling the parameters
 * given an allocation
 *
 * args:
 * - mu:    vector of location component
 * - m0:    mean of location's prior distribution (scalar)
 * - s20:   prior variance of the location component
 * - m1:    hyperparameter, mean of distribution of mu
 * - k1:    hyperparameter, scale factor of distribution of mu
 * - a1:    hyperparameter, shape parameter of distribution of s20
 * - b1:    hyperparameter, scale parameter of distribution of s20
 *
 * Void function
 ==================================================================================*/

void hyper_accelerate_MAR_L(arma::vec mu,
                            double &m0,
                            double &s20,
                            double m1,
                            double k1,
                            double a1,
                            double b1){


  double m_temp, k_temp, a_temp, b_temp, mu_m;
  int k = mu.n_elem;

  mu_m = arma::accu(mu) / k;
  k_temp = (k1 + k);
  m_temp = ((m1 * k1) + k * mu_m) / k_temp;
  a_temp = a1 + (k / 2.0);
  b_temp = b1 + (arma::accu(pow(mu - mu_m, 2)) + (k * k1 * pow(mu_m - m1, 2)) / (k_temp)) / 2;

  s20 = 1.0 / arma::randg(arma::distr_param(a_temp, 1 / b_temp));
  m0  = arma::randn() * sqrt(s20 / k_temp) + m_temp;

}

/*==================================================================================
 * Clean parameter - UNIVARIATE marginal sampler - LOCATION
 * discard the middle not used values for the clusters and
 * update the correspondent parameters.
 *
 * args:
 * - mu:      vector, each element a mean
 * - clust:   vector, each (integer) value is the cluster of the corresp. obs.
 *
 * Void function.
 ==================================================================================*/

void para_clean_MAR_L(arma::vec &mu,
                      arma::vec &clust){
  int k = mu.n_elem;
  double tmu;
  int u_bound;

  // for all the used parameters
  for(arma::uword i = 0; i < k; i++){

    // if a cluster is empty
    if((int) arma::sum(clust == i) == 0){

      // find the last full cluster, then swap
      for(arma::uword j = k; j > i; j--){
        if((int) arma::sum(clust == j) != 0){

          // SWAPPING!!
          clust( arma::find(clust == j) ).fill(i);

          tmu = mu[i];
          mu[i] = mu[j];
          mu[j] = tmu;

          break;
        }
      }
    }
  }

  // reduce dimensions
  u_bound = 0;
  for(arma::uword i = 0; i < k; i++){
    if(arma::accu(clust == i) > 0){
      u_bound += 1;
    }
  }

  // resize object to the correct dimension
  mu.resize(u_bound);
}

/*==================================================================================
 * Update clusters - UNIVARIATE marginal sampler - LOCATION SCALE
 *
 * args:
 * - data:        vector of observation
 * - mu:          mean values for each component
 * - s2:          variance
 * - clust:       vector of allocation
 * - mass:        mass paramter
 * - m0:          mean of location's prior distribution (scalar)
 * - s20:         prior variance of the location component
 * - a0:          shape of prior Gamma distribution on the scale component (scalar)
 * - b0:          rate of prior Gamma distribution on the scale component (scalar)
 * - iter:        current iteration
 * - new_val:     vector of new values
 * - sigma_PY:    discount parameter
 *
 * void function
 ==================================================================================*/

void clust_update_MAR_PY_L(arma::vec data,
                           arma::vec &mu,
                           double &s2,
                           arma::vec &clust,
                           double mass,
                           double m0,
                           double s20,
                           double a0,
                           double b0,
                           double sigma_PY){

  // initialize quantities
  int n = clust.n_elem;
  int k;
  arma::vec probs;
  double s_temp, m_temp;

  // loop over the observations
  for(arma::uword i = 0; i < n; i++){

    bool req_clean = false;
    if(arma::sum(clust == clust[i]) == 1){
      req_clean = true;
    }

    clust(i) = mu.n_elem + 1;
    if(req_clean){
      para_clean_MAR_L(mu,
                       clust);
    }

    // initialize useful quantities
    k = mu.n_elem;
    probs.resize(k+1);
    probs.fill(0);
    // arma::vec temp_clust = clust;
    // temp_clust(i) = k+1;

    // compute probabilities vector
    for(arma::uword j = 0; j < k; j++) {
      int nj = (int) arma::sum(clust == j);
      probs[j] = log(arma::normpdf(data[i], mu(j), sqrt(s2))) + log(nj - sigma_PY);
    }
    probs[k] = log(mass + k * sigma_PY) + log(arma::normpdf(data[i],
                                              (m0/s20 + data(i)/s2) * (1/ s20 + 1 / s2),
                                              1 / (1 / s20 + 1 / s2)));

    // sample new
    int temp_cl = rintnunif_log(probs);
    clust[i] = temp_cl;

    if(temp_cl == k){
      mu.resize(k+1);

      s_temp = 1 / (1/s2 + 1/s20);
      m_temp = s_temp * (m0 / s20 + data(i)/s2);

      mu(k) = arma::randn() * sqrt(s_temp) + m_temp;
    }
  }
}

/*----------------------------------------------------------------------
 *
 * UNIVARIATE
 * LOCATION-SCALE KERNEL
 * marginal functions
 *
 *----------------------------------------------------------------------
 */

/*==================================================================================
 * Accelerate - UNIVARIATE marginal sampler - LOCATION SCALE
 * acceleration step for reshuffling the parameters
 * given an allocation
 *
 * args:
 * - data:  vector of observation
 * - mu:    vector of location component
 * - s2:    vector of scale component
 * - clust: vector of allocation
 * - m0:    mean of distribution of location components
 * - k0:    scale factor of distribution of location components
 * - a0:    shape of distribution on the scale component (scalar)
 * - b0:    scale of distribution on the scale component (scalar)
 *
 * Void function
 ==================================================================================*/

void accelerate_MAR(arma::vec data,
                    arma::vec &mu,
                    arma::vec &s2,
                    arma::vec clust,
                    double m0,
                    double k0,
                    double a0,
                    double b0){
  int nj;
  double kn, mn, an, bn, data_m;
  arma::vec tdata;

  for(arma::uword j = 0; j < mu.n_elem; j++){
    nj = sum(clust == j);
    tdata = data.elem(arma::find(clust == j));

    data_m = arma::accu(data.elem(arma::find(clust == j))) / nj;
    kn = (k0 + nj);
    mn = ((m0 * k0) + nj * data_m) / kn;
    an = a0 + (nj / 2.0);
    bn = b0 + (arma::accu(pow(tdata - data_m, 2)) + (nj * k0 * pow(data_m - m0, 2)) / (kn)) / 2;

    s2[j] = 1.0 / arma::randg(arma::distr_param(an, 1 / bn));
    mu[j] = arma::randn() * sqrt(s2[j] / kn) + mn;
  }
}

/*==================================================================================
 * Hyper-accelerate - UNIVARIATE marginal sampler - LOCATION-SCALE
 * hyper-acceleration step for reshuffling the parameters
 * given an allocation
 *
 * args:
 * - mu:    vector of location component
 * - s2:    vector of scale component
 * - m0:    mean of distribution of location components
 * - k0:    scale factor of distribution of location components
 * - a0:    shape of distribution on the scale component
 * - b0:    scale of distribution on the scale component
 * - m1:    hyperparameter, location component of m0
 * - s21:   hyperparameter, scale component of m0
 * - tau1:  hyperparameter, shape component of k0
 * - tau2:  hyperparameter, rate component of k0
 * - a1:    hyperparameter, shape component of b0
 * - b1:    hyperparameter, rate component of b0
 *
 * Void function
 ==================================================================================*/

void hyper_accelerate_MAR(arma::vec mu,
                          arma::vec s2,
                          double &m0,
                          double &k0,
                          double a0,
                          double &b0,
                          double m1,
                          double s21,
                          double tau1,
                          double tau2,
                          double a1,
                          double b1){

  double m_temp, s2_temp, tau1_temp, tau2_temp, a_temp, b_temp, mu_m;
  int k = mu.n_elem;
  mu_m = arma::accu(mu) / k;

  tau1_temp = tau1 + k / 2;
  tau2_temp = tau2 + arma::accu(pow(mu - m0, 2) / s2) / 2;
  k0 = arma::randg(arma::distr_param(tau1_temp, 1 / tau2_temp));

  s2_temp = 1 / ( 1 / s21 + k0 * arma::accu(1 / s2) );
  m_temp  = s2_temp * ( m1 / s21 + k0 * arma::accu(mu / s2) );
  m0  = arma::randn() * sqrt(s2_temp) + m_temp;

  a_temp = a1 + k * a0;
  b_temp = b1 + arma::accu(1 / s2);
  b0 = arma::randg(arma::distr_param(a_temp, 1 / b_temp));
}

/*==================================================================================
 * Clean parameter - UNIVARIATE marginal sampler - LOCATION SCALE
 * discard the middle not used values for the clusters and
 * update the correspondent parameters.
 *
 * args:
 * - mu:      vector, each element a mean
 * - s2:      vector, each element a variance
 * - clust:   vector, each (integer) value is the cluster of the corresp. obs.
 *
 * Void function.
 ==================================================================================*/

void para_clean_MAR(arma::vec &mu,
                    arma::vec &s2,
                    arma::vec &clust) {
  int k = mu.n_elem;
  double tmu, ts2;
  int u_bound;

  // for all the used parameters
  for(arma::uword i = 0; i < k; i++){

    // if a cluster is empty
    if((int) arma::sum(clust == i) == 0){

      // find the last full cluster, then swap
      for(arma::uword j = k; j > i; j--){
        if((int) arma::sum(clust == j) != 0){

          // SWAPPING!!
          clust( arma::find(clust == j) ).fill(i);

          tmu = mu[i];
          mu[i] = mu[j];
          mu[j] = tmu;

          ts2 = s2[i];
          s2[i] = s2[j];
          s2[j] = ts2;

          break;
        }
      }
    }
  }

  // reduce dimensions
  u_bound = 0;
  for(arma::uword i = 0; i < k; i++){
    if(arma::accu(clust == i) > 0){
      u_bound += 1;
    }
  }

  // resize object to the correct dimension
  mu.resize(u_bound);
  s2.resize(u_bound);
}

/*==================================================================================
 * Update clusters - UNIVARIATE marginal sampler- LOCATION SCALE
 *
 * args:
 * - data:        vector of observation
 * - mu:          mean values for each component
 * - s2:          variance values for each component
 * - clust:       vector of allocation
 * - mass:        mass parameter
 * - m0:          mean of distribution of location components
 * - k0:          scale factor of distribution of location components
 * - a0:          shape of distribution on the scale component
 * - b0:          scale of distribution on the scale component
 * - iter:        current iteration
 * - new_val:     vector of new values
 * - sigma_PY:    discount parameter
 *
 * void function
 ==================================================================================*/

void clust_update_MAR_PY(arma::vec data,
                         arma::vec &mu,
                         arma::vec &s2,
                         arma::vec &clust,
                         double mass,
                         double m0,
                         double k0,
                         double a0,
                         double b0,
                         double sigma_PY){

  // initialize quantities
  int n = clust.n_elem;
  int k;
  arma::vec probs;
  double kn, mn, an, bn;

  // loop over the observations
  for(arma::uword i = 0; i < n; i++){

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
    k = mu.n_elem;
    probs.resize(k+1);
    probs.fill(0);
    // arma::vec temp_clust = clust;
    // temp_clust(i) = k+1;

    // compute probabilities vector
    for(arma::uword j = 0; j < k; j++) {
      int nj = (int) arma::sum(clust == j);
      probs[j] = log(arma::normpdf(data[i], mu(j), sqrt(s2(j)))) + log(nj - sigma_PY);
    }
    probs[k] = log(mass + k * sigma_PY) + dt_ls(data[i], 2 * a0, m0, sqrt(b0 * (1 + k0) / (a0*k0)));

    // sample new
    int temp_cl = rintnunif_log(probs);
    clust[i] = temp_cl;

    if(temp_cl == k){
      mu.resize(k+1);
      s2.resize(k+1);

      kn = (k0 + 1.0);
      mn = ((m0 * k0) + data[i]) / kn;
      an = a0 + (1.0 / 2.0);
      bn = b0 + ((k0 * pow(data(i) - m0, 2)) / (kn)) / 2;

      s2(k) = 1.0 / arma::randg(arma::distr_param(an, 1.0 / bn));
      mu(k) = arma::randn() * sqrt(s2(k) / kn) + mn;

    }
  }
}


/*----------------------------------------------------------------------
 *
 * MULTIVARIATE
 * LOCATION KERNEL
 * marginal functions
 *
 *----------------------------------------------------------------------
 */

/*==================================================================================
 * Accelerate - MULTIVARIATE marginal sampler - LOCATION
 * acceleration step for reshuffling the parameters
 * given an allocation
 *
 * args:
 * - data:  matrix of observations
 * - mu:    matrix of location components
 * - s2:    matrix - scale component
 * - clust: vector of allocations
 * - m0:    mean of location's prior distribution
 * - S20:   covariance matrix of the location component distribution
 * - S0:    matrix of Inverse Wishart distribution on s2
 * - n0:    degree of freedom of Inverse Wishart distribution on s2
 *
 * Void function
 ==================================================================================*/

void accelerate_MAR_mv_L(arma::mat data,
                         arma::mat &mu,
                         arma::mat &s2,
                         arma::vec clust,
                         arma::vec m0,
                         arma::mat S20,
                         arma::mat S0,
                         double n0){
  arma::mat S2_temp = s2;
  arma::vec m_temp = m0;
  arma::mat tdata, cdata, S_update(s2);
  S_update.fill(0);
  int nj;

  // loop over the different clusters
  for(arma::uword j = 0; j < mu.n_rows; j++){
    // initialize itra cluster quantities
    nj = sum(clust == j);
    tdata = data.rows(arma::find(clust == j));

    // update parameters
    S2_temp = arma::inv(arma::inv(S20) + nj * arma::inv(s2));
    m_temp = S2_temp * (arma::inv(s2) * arma::trans(sum(tdata, 0)) + arma::inv(S20) * m0);

    // sample from the posterior distributions
    mu.row(j)   = arma::trans(arma::mvnrnd(m_temp, S2_temp));

    cdata = tdata - arma::repmat(mu.row(j), nj, 1);
    S_update += arma::trans(cdata) * cdata;
  }

  int n = data.n_rows;
  arma::mat Sn = S0 + S_update;
  s2 = arma::inv(arma::wishrnd(arma::inv(Sn), (n0 + n)));
}

/*==================================================================================
 * Hyper-accelerate - MULTIVARIATE marginal sampler - LOCATION
 *
 * args:
 * - mu:      vector of location component
 * - m0:      mean of location's prior distribution (scalar)
 * - S20:     variance of the location component
 * - m1:      hyperparameter, location component of m0
 * - k1:      hyperparameter, scale parameter for variance of m0
 * - theta1:  hyperparameter, df of S20 distribution
 * - Theta2:  hyperparameter, matrix of S20 distribution
 *
 * Void function
 ==================================================================================*/

void hyper_accelerate_MAR_mv_L(arma::mat mu,
                               arma::vec &m0,
                               arma::mat &S20,
                               arma::vec m1,
                               double k1,
                               double theta1,
                               arma::mat Theta1){


  // compute the updated parameters
  arma::vec mum = arma::trans(mean(mu, 0));
  double k1n = k1 + mu.n_rows;
  arma::vec m1n = ((m1 * k1) + arma::trans(sum(mu, 0)))/k1n;
  arma::mat cmu = mu - arma::repmat(mum.t(), mu.n_rows, 1);
  double theta1n = theta1 + mu.n_rows;
  arma::mat Theta1n = Theta1 + arma::trans(cmu) * cmu + ((k1*mu.n_rows)/k1n) * (mum - m0) * arma::trans(mum - m0);

  // sample from the posterior distributions
  S20 = arma::inv(arma::wishrnd(arma::inv(Theta1n), theta1n));
  m0  = arma::mvnrnd(m1n, S20/k1n);
}

/*==================================================================================
 * Clean parameter - MULTIVARIATE marginal sampler - LOCATION
 * discard the middle not used values for the clusters and
 * update the correspondent parameters.
 *
 * args:
 * - mu:      matrix, each row a mean
 * - clust:   vector, each (integer) value is the cluster of the corresp. obs.
 *
 * Void function.
 ==================================================================================*/

void para_clean_MAR_mv_L(arma::mat &mu,
                         arma::vec &clust) {
  int k = mu.n_rows;
  int u_bound;

  // for all the used parameters
  for(arma::uword i = 0; i < k; i++){

    // if a cluster is empty
    if((int) arma::sum(clust == i) == 0){

      // find the last full cluster, then swap
      for(arma::uword j = k; j > i; j--){
        if((int) arma::sum(clust == j) != 0){

          // SWAPPING!!
          clust( arma::find(clust == j) ).fill(i);
          mu.swap_rows(i,j);
          break;
        }
      }
    }
  }

  // reduce dimensions
  u_bound = 0;
  for(arma::uword i = 0; i < k; i++){
    if(arma::accu(clust == i) > 0){
      u_bound += 1;
    }
  }

  // resize object to the correct dimension
  mu.resize(u_bound,mu.n_cols);
}

/*==================================================================================
 * Update clusters - MULTIVARIATE marginal sampler - LOCATION
 *
 * args:
 * - data:        vector of observation
 * - mu:          mean values for each component
 * - s2:          covariance
 * - clust:       vector of allocation
 * - mass:        mass parameter
 * - m0:          mean of location's prior distribution
 * - S20:         covariance matrix of the location component distribution
 * - S0:          matrix of Inverse Wishart distribution on s2
 * - n0:          degree of freedom of Inverse Wishart distribution on s2
 * - iter:        current iteration
 * - new_val:     vector of new values
 * - sigma_PY:    discount parameter
 *
 * void function
 ==================================================================================*/

void clust_update_MAR_PY_mv_L(arma::mat data,
                              arma::mat &mu,
                              arma::mat s2,
                              arma::vec &clust,
                              double mass,
                              arma::vec m0,
                              arma::mat S20,
                              arma::mat S0,
                              double n0,
                              double sigma_PY){
  // initialize quantities
  int n = clust.n_elem;
  int d = data.n_cols;
  int k;
  arma::vec probs, m_temp;
  arma::mat Sn, S2_temp;

  arma::mat rooti  = arma::trans(arma::inv(trimatu(arma::chol(s2))));
  double rsum = arma::sum(log(rooti.diag()));

  // loop over the observations
  for(arma::uword i = 0; i < n; i++){

    bool req_clean = false;
    if(arma::sum(clust == clust[i]) == 1){
      req_clean = true;
    }

    clust(i) = mu.n_rows + 1;
    if(req_clean){
      para_clean_MAR_mv_L(mu,
                          clust);
    }

    // initialize useful quantities
    k = mu.n_rows;
    probs.resize(k+1);
    probs.fill(0);

    // compute probabilities vector

    for(arma::uword j = 0; j < k; j++) {
      int nj = (int) arma::sum(clust == j);
      arma::vec cdata  = rooti * arma::trans(data.row(i) - mu.row(j)) ;
      double out       = - (d / 2.0) * log(2.0 * M_PI) - 0.5 * arma::sum(cdata%cdata) + rsum;
      probs[j] = out + log(nj - sigma_PY);
    }

    S2_temp = arma::inv(arma::inv(s2) + arma::inv(S20));
    m_temp  = S2_temp * (arma::inv(s2) * arma::trans(data.row(i)) + arma::inv(S20) * m0);
    arma::mat rooti2  = arma::trans(arma::inv(trimatu(arma::chol(s2 + S2_temp))));
    arma::vec cdata2  = rooti2 * arma::trans(data.row(i) - m_temp.t()) ;
    double out2       = - (d / 2.0) * log(2.0 * M_PI) - 0.5 * arma::sum(cdata2%cdata2) + arma::sum(log(rooti2.diag()));
    probs[k] = log(mass + k * sigma_PY) + out2;

    // sample new
    int temp_cl = rintnunif_log(probs);
    clust[i] = temp_cl;

    if(temp_cl == k){
      mu.resize(k+1, mu.n_cols);
      mu.row(k)   = arma::trans(arma::mvnrnd(m_temp, S2_temp));
    }
  }
}

/*----------------------------------------------------------------------
 *
 * MULTIVARIATE
 * LOCATION-SCALE KERNEL
 * marginal functions
 *
 *----------------------------------------------------------------------
 */

/*==================================================================================
 * Accelerate - MULTIVARIATE marginal sampler - LOCATION SCALE
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
 ==================================================================================*/

void accelerate_MAR_mv(arma::mat data,
                       arma::mat &mu,
                       arma::cube &s2,
                       arma::vec clust,
                       arma::vec m0,
                       double k0,
                       arma::mat S0,
                       double n0){
  int nj;
  double kn, nn;
  arma::vec mn;
  arma::mat cdata;
  arma::mat Sn;

  // loop over the different clusters
  for(arma::uword j = 0; j < mu.n_rows; j++){
    // initialize itra cluster quantities
    nj = sum(clust == j);
    arma::mat tdata = data.rows(arma::find(clust == j));

    // update parameters
    kn = k0 + nj;
    mn = ((m0 * k0) + arma::trans(sum(tdata, 0)))/kn;
    cdata = tdata - arma::repmat(mean(tdata, 0), nj, 1);
    nn = n0 + nj;
    Sn = S0 + arma::trans(cdata) * cdata + ((k0*nj)/kn) *
      (arma::trans(mean(tdata, 0)) - m0) * arma::trans(arma::trans(mean(tdata, 0)) - m0);

    // sample from the posterior distributions
    s2.slice(j) = arma::inv(arma::wishrnd(arma::inv(Sn), nn));
    mu.row(j)   = arma::trans(arma::mvnrnd(mn, s2.slice(j)/kn));
  }
}

/*==================================================================================
 * Hyper-accelerate - MULTIVARIATE marginal sampler - LOCATION-SCALE
 *
 * args:
 * - mu:      vector of location component
 * - s2:      covariance matricies
 * - m0:      mean of location distribution
 * - k0:      scale factor of location distribution
 * - S20:     matrix of the scale distribution
 * - n0:      df of scale distribution
 * - m1:      hyperparameter, location component of m0
 * - k1:      hyperparameter, scale parameter for variance of m0
 * - tau1:    hyperparameter, shape of k0 distribution
 * - tau2:    hyperparameter, rate of k0 distribution
 * - theta1:  hyperparameter, df of S20 distribution
 * - Theta2:  hyperparameter, matrix of S20 distribution
 *
 * Void function
 ==================================================================================*/

void hyper_accelerate_MAR_mv_LS(arma::mat mu,
                                arma::cube s2,
                                arma::vec &m0,
                                double &k0,
                                arma::mat &S0,
                                double n0,
                                arma::vec m1,
                                arma::mat S1,
                                double tau1,
                                double tau2,
                                double theta1,
                                arma::mat Theta1){

  arma::mat accu_Sinv(S0.n_cols, S0.n_cols, arma::fill::zeros);
  arma::vec accu_Smu(m0.n_elem, arma::fill::zeros);
  double accu_qform = 0.0;

  for(arma::uword j = 0; j < mu.n_rows; j++){
    accu_Sinv   += arma::inv(s2.slice(j));
    accu_Smu    += arma::inv(s2.slice(j)) * mu.row(j).t();
    accu_qform  += arma::as_scalar(mu.row(j) * arma::inv(s2.slice(j)) * mu.row(j).t());
  }

  // update m0
  arma::mat S1n = arma::inv(arma::inv(S1) + k0 * accu_Sinv);
  arma::vec m1n = S1n * (arma::inv(S1) * m1 + k0 * accu_Smu);
  m0 = arma::mvnrnd(m1n, S1n);

  // update S0
  double theta1n = theta1 + mu.n_rows * n0;
  arma::mat Theta1n = arma::inv(arma::inv(Theta1) + accu_Sinv);
  S0 = arma::wishrnd(Theta1n, theta1n);

  // update k0
  double tau1n = tau1 + mu.n_rows;
  double tau2n = tau2 + accu_qform;
  k0 = arma::randg(arma::distr_param(tau1n, 1 / tau2n));
}

/*==================================================================================
 * Clean parameter - MULTIVARIATE marginal sampler - LOCATION SCALE
 * discard the middle not used values for the clusters and
 * update the correspondent parameters.
 *
 * args:
 * - mu:      matrix, each row a mean
 * - s2:      cube, each slice a covariance matrix
 * - clust:   vector, each (integer) value is the cluster of the corresp. obs.
 *
 * Void function.
 ==================================================================================*/

void para_clean_MAR_mv(arma::mat &mu,
                       arma::cube &s2,
                       arma::vec &clust) {
  int k = mu.n_rows;
  int u_bound;

  // for all the used parameters
  for(arma::uword i = 0; i < k; i++){

    // if a cluster is empty
    if((int) arma::sum(clust == i) == 0){

      // find the last full cluster, then swap
      for(arma::uword j = k; j > i; j--){
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
  u_bound = 0;
  for(arma::uword i = 0; i < k; i++){
    if(arma::accu(clust == i) > 0){
      u_bound += 1;
    }
  }

  // resize object to the correct dimension
  mu.resize(u_bound,mu.n_cols);
  s2.resize(s2.n_rows, s2.n_cols, u_bound);
}

/*==================================================================================
 * Update clusters - MULTIVARIATE importance conditional sampler - LOCATION SCALE
 *
 * args:
 * - data:        vector of observation
 * - mu:          mean values for each component
 * - s2:          covariance for each component
 * - clust:       vector of allocation
 * - mass:        mass parameter
 * - m0:          vector of location's prior distribution
 * - k0:          double of NIG on scale of location distribution
 * - S0:          matrix of Inverse Wishart distribution
 * - n0:          degree of freedom of Inverse Wishart distribution
 * - iter:        current iteration
 * - new_val:     vector of new values
 * - sigma_PY:    discount parameter
 *
 * void function
 ==================================================================================*/

void clust_update_MAR_PY_mv(arma::mat data,
                            arma::mat &mu,
                            arma::cube &s2,
                            arma::vec &clust,
                            double mass,
                            arma::vec m0,
                            double k0,
                            arma::mat S0,
                            double n0,
                            double sigma_PY){
  // initialize quantities
  int n = clust.n_elem;
  int d = data.n_cols;
  int k;
  arma::vec probs, mn;
  double kn, nn;
  arma::mat Sn;

  // loop over the observations
  for(arma::uword i = 0; i < n; i++){

    bool req_clean = (arma::sum(clust == clust[i]) == 1);

    clust(i) = mu.n_rows + 1;
    if(req_clean){
      para_clean_MAR_mv(mu,
                        s2,
                        clust);
    }

    // initialize useful quantities
    k = mu.n_rows;
    probs.resize(k+1);
    probs.fill(0);

    // compute probabilities vector
    for(arma::uword j = 0; j < k; j++) {
      int nj = (int) arma::sum(clust == j);
      arma::mat rooti  = arma::trans(arma::inv(trimatu(arma::chol(s2.slice(j)))));
      arma::vec cdata  = rooti * arma::trans(data.row(i) - mu.row(j)) ;
      double out       = - (d / 2.0) * log(2.0 * M_PI) - 0.5 * arma::sum(cdata%cdata) + arma::sum(log(rooti.diag()));
      probs(j) = out + log(nj - sigma_PY);
    }
    probs(k) = log(mass + k * sigma_PY) + dt_ls_mv(data.row(i).t(), n0 - d + 1, m0, S0 * (k0 + 1) / (k0 * (n0 - d + 1)));

    // sample new
    int temp_cl = rintnunif_log(probs);
    clust(i) = temp_cl;

    if(temp_cl == k){
      mu.resize(k+1, mu.n_cols);
      s2.resize(s2.n_rows, s2.n_cols, k+1);

      // update parameters
      kn = k0 + 1;
      mn = ((m0 * k0) + data.row(i).t())/kn;
      nn = n0 + 1;

      Sn = S0 + (data.row(i).t() - mn) * arma::trans(data.row(i).t() - mn) + ((k0+1)/kn) *
        (data.row(i).t() - m0) * arma::trans(data.row(i).t() - m0);

      // sample from the posterior distributions
      s2.slice(k) = arma::inv(arma::wishrnd(arma::inv(Sn), nn));
      mu.row(k)   = arma::trans(arma::mvnrnd(mn, s2.slice(k)/kn));
    }
  }
}

/*----------------------------------------------------------------------
 *
 * MULTIVARIATE
 * PRODUCT KERNEL
 * marginal functions
 *
 *----------------------------------------------------------------------
 */

/*==================================================================================
 * Accelerate - MULTIVARIATE marginal sampler - PRODUCT
 * acceleration step for reshuffling the parameters
 * given an allocation
 *
 * args:
 * - data:  matrix of observation
 * - mu:    matrix of location component
 * - s2:    matrix of scale component
 * - clust: vector of allocation
 * - m0:    vector of location's prior distribution
 * - k0:    vector of location's variance tuning parameter, one for each dimension
 * - a0:    vector of shape parameters for the scale component, one for each dimension
 * - b0:    vector of scale parameters for the scale component, one for each dimension
 *
 * Void function
 ==================================================================================*/

void accelerate_MAR_mv_P(arma::mat data,
                         arma::mat &mu,
                         arma::mat &s2,
                         arma::vec clust,
                         arma::vec m0,
                         arma::vec k0,
                         arma::vec a0,
                         arma::vec b0){
  double an, bn, kn, mn, data_m;
  arma::mat cdata;
  int nj;

  // loop over the different clusters
  for(arma::uword j = 0; j < mu.n_rows; j++){

    nj = sum(clust == j);
    arma::mat tdata = data.rows(arma::find(clust == j));

    for(arma::uword l = 0; l < mu.n_cols; l++){

      data_m = arma::accu(tdata.col(l)) / nj;
      kn = k0(l) + nj;
      mn = ((m0(l) * k0(l)) + nj * data_m) / kn;
      an = a0(l) + (nj / 2.0);
      bn = b0(l) + (arma::accu(pow(tdata.col(l) - data_m, 2)) + (nj * k0(l) * pow(data_m - m0(l), 2)) / (kn)) / 2;

      s2(j,l) = 1.0 / arma::randg(arma::distr_param(an, 1.0 / bn));
      mu(j,l) = arma::randn() * sqrt(s2(j,l) / kn) + mn;
    }
  }
}

/*==================================================================================
 * Hyper-accelerate - MULTIVARIATE marginal sampler - PRODUCT
 * acceleration step for reshuffling the parameters
 * given an allocation
 *
 * args:
 * - mu:    vector of location component
 * - s2:    vector of scale component
 * - m0:    mean of location's prior distribution (scalar)
 * - k0:    tuning parameter of variance of the location component
 * - a0:    shape parameters of scale component
 * - b0:    scale parameter of scale component
 * - m1:    hyperparameter, location component of m0
 * - s21:   hyperparameter, scale component of m0
 * - tau1:  hyperparameter, shape component of k0
 * - tau2:  hyperparameter, rate component of k0
 * - a1:    hyperparameter, shape component of b0
 * - b1:    hyperparameter, rate component of b0
 *
 * Void function
 ==================================================================================*/

void hyper_accelerate_MAR_mv_P(arma::mat mu,
                               arma::mat s2,
                               arma::vec &m0,
                               arma::vec &k0,
                               arma::vec a0,
                               arma::vec &b0,
                               arma::vec m1,
                               arma::vec s21,
                               arma::vec tau1,
                               arma::vec tau2,
                               arma::vec a1,
                               arma::vec b1){

  double m_temp, s2_temp, tau1_temp, tau2_temp, a_temp, b_temp, mu_m;
  int k = mu.n_rows;

  for(arma::uword j = 0; j < mu.n_cols; j++){
    mu_m = arma::accu(mu.col(j)) / k;

    tau1_temp = tau1(j) + k / 2;
    tau2_temp = tau2(j) + arma::accu(pow(mu.col(j) - m0(j), 2) / s2.col(j)) / 2;
    k0(j) = arma::randg(arma::distr_param(tau1_temp, 1 / tau2_temp));

    s2_temp = 1 / ( 1 / s21(j) + k0(j) * arma::accu(1 / s2.col(j)) );
    m_temp  = s2_temp * ( m1(j) / s21(j) + k0(j) * arma::accu(mu.col(j) / s2.col(j)) );
    m0(j)  = arma::randn() * sqrt(s2_temp) + m_temp;

    a_temp = a1(j) + k * a0(j);
    b_temp = b1(j) + arma::accu(1 / s2.col(j));
    b0(j) = arma::randg(arma::distr_param(a_temp, 1 / b_temp));
  }
}


/*==================================================================================
 * Clean parameter - MULTIVARIATE marginal sampler - PRODUCT
 * discard the middle not used values for the clusters and
 * update the correspondent parameters.
 *
 * args:
 * - mu:      matrix, each row a mean
 * - s2:      matrix, each row a variance vector
 * - clust:   vector, each (integer) value is the cluster of the corresp. obs.
 *
 * Void function.
 ==================================================================================*/

void para_clean_MAR_mv_P(arma::mat &mu,
                         arma::mat &s2,
                         arma::vec &clust) {
  int k = mu.n_rows;
  int u_bound;

  // for all the used parameters
  for(arma::uword i = 0; i < k; i++){

    // if a cluster is empty
    if((int) arma::sum(clust == i) == 0){

      // find the last full cluster, then swap
      for(arma::uword j = k; j > i; j--){
        if((int) arma::sum(clust == j) != 0){

          // SWAPPING!!
          clust( arma::find(clust == j) ).fill(i);
          mu.swap_rows(i,j);
          s2.swap_rows(i,j);
          break;
        }
      }
    }
  }

  // reduce dimensions
  u_bound = 0;
  for(arma::uword i = 0; i < k; i++){
    if(arma::accu(clust == i) > 0){
      u_bound += 1;
    }
  }

  // resize object to the correct dimension
  mu.resize(u_bound, mu.n_cols);
  s2.resize(u_bound, s2.n_cols);
}

/*==================================================================================
 * Update clusters - MULTIVARIATE marginal sampler - PRODUCT
 *
 * args:
 * - data:        vector of observation
 * - mu:          mean values for each component
 * - s2:          variance values for each component
 * - clust:       vector of allocation
 * - mass:        mass parameter
 * - m0:          vector of location's prior distribution
 * - k0:          vector of location's variance tuning parameter, one for each dimension
 * - a0:          vector of shape parameters for the scale component, one for each dimension
 * - b0:          vector of scale parameters for the scale component, one for each dimension
 * - iter:        current iteration
 * - new_val:     vector of new values
 * - sigma_PY:    discount parameter
 *
 * void function
 ==================================================================================*/

void clust_update_MAR_PY_mv_P(arma::mat data,
                              arma::mat &mu,
                              arma::mat &s2,
                              arma::vec &clust,
                              double mass,
                              arma::vec m0,
                              arma::vec k0,
                              arma::vec a0,
                              arma::vec b0,
                              double sigma_PY){
  // initialize quantities
  int n = clust.n_elem;
  int d = data.n_cols;
  int k;
  arma::vec probs;
  double kn, an, bn, mn;

  // loop over the observations
  for(arma::uword i = 0; i < n; i++){

    bool req_clean = false;
    if(arma::sum(clust == clust[i]) == 1){
      req_clean = true;
    }

    clust(i) = mu.n_rows + 1;
    if(req_clean){
      para_clean_MAR_mv_P(mu,
                          s2,
                          clust);
    }

    // initialize useful quantities
    k = mu.n_rows;
    probs.resize(k+1);
    probs.fill(0);

    // compute probabilities vector
    for(arma::uword j = 0; j < k; j++) {

      probs(j) = log(arma::sum(clust == j)- sigma_PY);
      for(arma::uword l = 0; l < d; l++){
        probs(j) +=  - log(2 * M_PI * s2(j,l)) / 2 - (pow(data(i,l) - mu(j,l), 2) / (2 * s2(j,l)));
      }
    }
    probs(k) = log(mass + k * sigma_PY);
    for(arma::uword l = 0; l < d; l++){
      probs(k) +=  dt_ls(data(i,l), 2 * a0(l), m0(l), sqrt(b0(l) * (1 + k0(l)) / (k0(l) * a0(l))));
    }

    // sample new
    int temp_cl = rintnunif_log(probs);
    clust[i] = temp_cl;

    if(temp_cl == k){
      mu.resize(k+1, mu.n_cols);
      s2.resize(k+1, s2.n_cols);

      for(arma::uword l = 0; l < d; l++){
        kn = k0(l) + 1;
        mn = ((m0(l) * k0(l)) + data(i,l)) / kn;
        an = a0(l) + (1 / 2.0);
        bn = b0(l) + ((k0(l) * pow(data(i,l) - m0(l), 2)) / (kn)) / 2;

        s2(k,l) = 1.0 / arma::randg(arma::distr_param(an, 1.0 / bn));
        mu(k,l) = arma::randn() * sqrt(s2(k,l) / kn) + mn;
      }
    }
  }
}

/*----------------------------------------------------------------------
 *
 * MIXTURE OF REGRESSION KERNELS
 * LOCATION-SCALE KERNEL
 * MAR functions
 *
 *----------------------------------------------------------------------
 */

/*==================================================================================
 * Accelerate - MULTIVARIATE importance conditional sampler - MRK
 * acceleration step for reshuffling the parameters
 * given an allocation
 *
 * args:
 * - y:       target variable (n x 1)
 * - covs:    covariates (n x (d + 1))
 * - beta:    regression coefficients (k x (d+1))
 * - sigma2:  scale component of kernel function
 * - clust:   vector of allocation
 * - beta0:   vector of location's prior distribution
 * - Sb0:     double of NIG on scale of location distribution
 * - a0:      shape parameter of scale component
 * - b0:      scale parameter of scale component
 *
 * Void function
 ==================================================================================*/

void accelerate_MAR_mv_MRK(arma::vec y,
                           arma::mat covs,
                           arma::mat &beta,
                           arma::vec &sigma2,
                           arma::vec clust,
                           arma::vec beta0,
                           arma::mat Sb0,
                           double a0,
                           double b0){
  arma::mat cdata;
  arma::mat tdata;
  arma::vec ty;

  double an, bn;
  arma::mat tSb;
  arma::vec tbeta0;

  // loop over the different clusters
  for(arma::uword j = 0; j < beta.n_rows; j++){

    // initialize itra cluster quantities
    int nj = sum(clust == j);
    tdata = covs.rows(arma::find(clust == j));
    ty = y.elem(arma::find(clust == j));
    cdata = ty - (tdata * beta.row(j).t());

    // update parameters
    tSb = arma::inv(arma::inv(Sb0) + arma::trans(tdata) * tdata / sigma2(j));
    tbeta0 = tSb * (arma::inv(Sb0) * beta0 + (arma::trans(tdata) * ty) / sigma2(j));
    an = a0 + nj/2;
    bn = b0 + (arma::accu(pow(cdata, 2))) / 2;

    // sample from the posterior distributions
    sigma2(j) = 1.0 / arma::randg(arma::distr_param(an, 1.0 / (bn)));
    beta.row(j) = arma::trans(arma::mvnrnd(tbeta0, tSb));
  }
}

/*==================================================================================
 * Hyper-accelerate - MULTIVARIATE importance conditional sampler - MRK
 *
 * args:
 * - mu:      vector of location component
 * - s2:      covariance matricies
 * - m0:      mean of location distribution
 * - k0:      scale factor of location distribution
 * - S20:     matrix of the scale distribution
 * - n0:      df of scale distribution
 * - m1:      hyperparameter, location component of m0
 * - k1:      hyperparameter, scale parameter for variance of m0
 * - tau1:    hyperparameter, shape of k0 distribution
 * - tau2:    hyperparameter, rate of k0 distribution
 * - theta1:  hyperparameter, df of S20 distribution
 * - Theta2:  hyperparameter, matrix of S20 distribution
 *
 * Void function
 ==================================================================================*/

void hyper_accelerate_MAR_mv_MRK(arma::vec y,
                                 arma::mat covs,
                                 arma::vec clust,
                                 arma::mat beta,
                                 arma::vec sigma2,
                                 arma::vec &beta0,
                                 arma::mat &Sb0,
                                 double a0,
                                 double &b0,
                                 arma::vec beta1,
                                 double k1,
                                 double sb1,
                                 arma::mat Sb1,
                                 double tau1,
                                 double tau2){

  // initialize quantities
  arma::mat cbeta = beta - arma::repmat(mean(beta, 0), beta.n_rows, 1);

  // sampling hyperparameters
  double kn = k1 + beta.n_rows;
  arma::vec betan = ((beta1 * k1) + arma::trans(sum(beta, 0)))/kn;
  double sbn = sb1 + beta.n_rows;
  arma::mat Sbn = Sb1 + arma::trans(cbeta) * cbeta + ((k1 * beta.n_rows) / kn) *
    (arma::trans(mean(beta, 0)) - beta1) * arma::trans(arma::trans(mean(beta, 0)) - beta1);
  double tau1n = tau1 + beta.n_rows * a0;
  double tau2n = tau2 + arma::accu(1 / sigma2);

  Sb0 = arma::inv(arma::wishrnd(arma::inv(Sbn), sbn));
  beta0 = arma::mvnrnd(betan, Sb0 / kn);
  b0 = arma::randg(arma::distr_param(tau1n, 1 / tau2n));

}


/*==================================================================================
 * Clean parameter - MULTIVARIATE importance conditional sampler - MRK
 * discard the middle not used values for the clusters and
 * update the correspondent parameters.
 *
 * args:
 * - mu:      matrix, each row a mean
 * - s2:      cube, each slice a covariance matrix
 * - clust:   vector, each (integer) value is the cluster of the corresp. obs.
 *
 * Void function.
 ==================================================================================*/

void para_clean_MAR_mv_MRK(arma::mat &beta,
                           arma::vec &sigma2,
                           arma::vec &clust) {
  int k = beta.n_rows;
  double tsigma2;

  // for all the used parameters
  for(arma::uword i = 0; i < k; i++){

    // if a cluster is empty
    if((int) arma::sum(clust == i) == 0){

      // find the last full cluster, then swap
      for(arma::uword j = k; j > i; j--){
        if((int) arma::sum(clust == j) != 0){

          // swap the corresponding elements
          clust( arma::find(clust == j) ).fill(i);
          beta.swap_rows(i,j);
          tsigma2 = sigma2(i);
          sigma2(i) = sigma2(j);
          sigma2(j) = tsigma2;
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
  beta.resize(u_bound, beta.n_cols);
  sigma2.resize(u_bound);
}

/*==================================================================================
 * Update clusters - MULTIVARIATE importance conditional sampler - MRK
 *
 * args:
 * - data:        vector of observation
 * - mujoin:      mean values for each component
 * - s2join:      variance for each component
 * - probjoin:    prob for each component
 * - clust:       vector of allocation
 * - max_val:     vector, number of already existent atoms
 * - iter:        current iteration
 * - new_val:     vector of new values
 *
 * void function
 ==================================================================================*/

void clust_update_MAR_mv_MRK(arma::vec y,
                             arma::mat covs,
                             arma::mat &beta,
                             arma::vec &sigma2,
                             arma::vec &clust,
                             double mass,
                             arma::vec beta0,
                             arma::mat Sb0,
                             double a0,
                             double b0,
                             double sigma_PY,
                             int napprox){

  // initialize quantities
  int n = clust.n_elem;
  int k;
  arma::vec probs, betan, sigma2temp;
  double an, bn, prob_temp;
  arma::mat Sbn, betatemp;

  // sample the temp values
  sigma2temp = 1.0 / arma::randg(napprox, arma::distr_param(a0, 1.0 / b0));
  betatemp = arma::trans(arma::mvnrnd(beta0, Sb0, napprox));

  // loop over the observations
  for(arma::uword i = 0; i < n; i++){

    bool req_clean = false;
    if(arma::sum(clust == clust[i]) == 1){
      req_clean = true;
    }

    clust(i) = beta.n_rows + 1;
    if(req_clean){
      para_clean_MAR_mv_MRK(beta,
                            sigma2,
                            clust);
    }

    // initialize useful quantities
    k = beta.n_rows;
    probs.resize(k+1);
    probs.fill(0);

    // compute probabilities vector
    for(arma::uword j = 0; j < k; j++) {
      int nj = (int) arma::sum(clust == j);
      probs(j) = log(nj - sigma_PY) - log(2 * M_PI * sigma2(j)) / 2 -
        (pow(y(i) - arma::dot(covs.row(i), beta.row(j)), 2) / (2 * sigma2(j)));
    }

    prob_temp = 0.0;
    for(arma::uword l = 0; l < napprox; l++){
      prob_temp += arma::normpdf(y(i), arma::dot(covs.row(i), betatemp.row(l)), sqrt(sigma2temp(l)));
    }
    probs(k) = log(mass + k * sigma_PY) + log(prob_temp / napprox);

    // sample new
    int temp_cl = rintnunif_log(probs);
    clust(i) = temp_cl;

    if(temp_cl == k){

      beta.resize(k+1, beta.n_cols);
      sigma2.resize(k+1);

      an = a0 + (1.0 / 2.0);
      bn = b0 + pow(y(i) - arma::dot(covs.row(i), beta0), 2) / 2;
      sigma2(k) = 1.0 / arma::randg(arma::distr_param(an, 1.0 / (bn)));

      Sbn = arma::inv(arma::inv(Sb0) + arma::trans(covs.row(i)) * covs.row(i) / sigma2(k));
      betan = Sbn * (arma::inv(Sb0) * beta0 + (arma::trans(covs.row(i)) * y(i)) / sigma2(k));
      beta.row(k) = arma::trans(arma::mvnrnd(betan, Sbn));

    }
  }
}

/*----------------------------------------------------------------------
 *
 * MIXTURE OF REGRESSION KERNELS - LOCATION
 * LOCATION-SCALE KERNEL
 * MAR functions
 *
 *----------------------------------------------------------------------
 */

/*==================================================================================
 * Accelerate - MULTIVARIATE importance conditional sampler - MRK
 * acceleration step for reshuffling the parameters
 * given an allocation
 *
 * args:
 * - y:       target variable (n x 1)
 * - covs:    covariates (n x (d + 1))
 * - beta:    regression coefficients (k x (d+1))
 * - sigma2:  scale component of kernel function
 * - clust:   vector of allocation
 * - beta0:   vector of location's prior distribution
 * - Sb0:     double of NIG on scale of location distribution
 * - a0:      shape parameter of scale component
 * - b0:      scale parameter of scale component
 *
 * Void function
 ==================================================================================*/

void accelerate_MAR_mv_MRK_L(arma::vec y,
                             arma::mat covs,
                             arma::mat &beta,
                             double &sigma2,
                             arma::vec clust,
                             arma::vec beta0,
                             arma::mat Sb0,
                             double a0,
                             double b0){
  arma::mat cdata;
  arma::mat tdata;
  arma::vec ty;

  double an, bn;
  arma::mat tSb;
  arma::vec tbeta0;

  double accu_cdata = 0.0;

  // loop over the different clusters
  for(arma::uword j = 0; j < beta.n_rows; j++){

    // initialize itra cluster quantities
    int nj = sum(clust == j);
    tdata = covs.rows(arma::find(clust == j));
    ty = y.elem(arma::find(clust == j));
    cdata = ty - (tdata * beta.row(j).t());

    accu_cdata += arma::accu(pow(cdata, 2));

    // update parameters
    tSb = arma::inv(arma::inv(Sb0) + arma::trans(tdata) * tdata / sigma2);
    tbeta0 = tSb * (arma::inv(Sb0) * beta0 + (arma::trans(tdata) * ty) / sigma2);
    beta.row(j) = arma::trans(arma::mvnrnd(tbeta0, tSb));
  }

  // sample from the posterior distributions
  an = a0 + y.n_elem / 2;
  bn = b0 + accu_cdata / 2;
  sigma2 = 1.0 / arma::randg(arma::distr_param(an, 1.0 / (bn)));

}

/*==================================================================================
 * Hyper-accelerate - MULTIVARIATE importance conditional sampler - MRK
 *
 * args:
 * - mu:      vector of location component
 * - s2:      covariance matricies
 * - m0:      mean of location distribution
 * - k0:      scale factor of location distribution
 * - S20:     matrix of the scale distribution
 * - n0:      df of scale distribution
 * - m1:      hyperparameter, location component of m0
 * - k1:      hyperparameter, scale parameter for variance of m0
 * - theta1:  hyperparameter, df of S20 distribution
 * - Theta2:  hyperparameter, matrix of S20 distribution
 *
 * Void function
 ==================================================================================*/

void hyper_accelerate_MAR_mv_MRK_L(arma::vec y,
                                   arma::mat covs,
                                   arma::vec clust,
                                   arma::mat beta,
                                   arma::vec &beta0,
                                   arma::mat &Sb0,
                                   double a0,
                                   double &b0,
                                   arma::vec beta1,
                                   double k1,
                                   double sb1,
                                   arma::mat Sb1){

  // initialize quantities
  arma::mat cbeta = beta - arma::repmat(mean(beta, 0), beta.n_rows, 1);

  // sampling hyperparameters
  double kn = k1 + beta.n_rows;
  arma::vec betan = ((beta1 * k1) + arma::trans(sum(beta, 0)))/kn;
  double sbn = sb1 + beta.n_rows;
  arma::mat Sbn = Sb1 + arma::trans(cbeta) * cbeta + ((k1 * beta.n_rows) / kn) *
    (arma::trans(mean(beta, 0)) - beta1) * arma::trans(arma::trans(mean(beta, 0)) - beta1);

  Sb0 = arma::inv(arma::wishrnd(arma::inv(Sbn), sbn));
  beta0 = arma::mvnrnd(betan, Sb0 / kn);

}


/*==================================================================================
 * Clean parameter - MULTIVARIATE importance conditional sampler - MRK
 * discard the middle not used values for the clusters and
 * update the correspondent parameters.
 *
 * args:
 * - mu:      matrix, each row a mean
 * - s2:      cube, each slice a covariance matrix
 * - clust:   vector, each (integer) value is the cluster of the corresp. obs.
 *
 * Void function.
 ==================================================================================*/

void para_clean_MAR_mv_MRK_L(arma::mat &beta,
                             arma::vec &clust) {
  int k = beta.n_rows;
  double tsigma2;

  // for all the used parameters
  for(arma::uword i = 0; i < k; i++){

    // if a cluster is empty
    if((int) arma::sum(clust == i) == 0){

      // find the last full cluster, then swap
      for(arma::uword j = k; j > i; j--){
        if((int) arma::sum(clust == j) != 0){

          // swap the corresponding elements
          clust( arma::find(clust == j) ).fill(i);
          beta.swap_rows(i,j);
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
  beta.resize(u_bound, beta.n_cols);
}

/*==================================================================================
 * Update clusters - MULTIVARIATE importance conditional sampler - MRK
 *
 * args:
 * - data:        vector of observation
 * - mujoin:      mean values for each component
 * - s2join:      variance for each component
 * - probjoin:    prob for each component
 * - clust:       vector of allocation
 * - max_val:     vector, number of already existent atoms
 * - iter:        current iteration
 * - new_val:     vector of new values
 *
 * void function
 ==================================================================================*/

void clust_update_MAR_mv_MRK_L(arma::vec y,
                               arma::mat covs,
                               arma::mat &beta,
                               double sigma2,
                               arma::vec &clust,
                               double mass,
                               arma::vec beta0,
                               arma::mat Sb0,
                               double sigma_PY,
                               int napprox){

  // initialize quantities
  int n = clust.n_elem;
  int k;
  arma::vec probs, betan;
  double an, bn, prob_temp;
  arma::mat Sbn, betatemp;

  // sample the temp values
  betatemp = arma::trans(arma::mvnrnd(beta0, Sb0, napprox));

  // loop over the observations
  for(arma::uword i = 0; i < n; i++){

    bool req_clean = false;
    if(arma::sum(clust == clust[i]) == 1){
      req_clean = true;
    }

    clust(i) = beta.n_rows + 1;
    if(req_clean){
      para_clean_MAR_mv_MRK_L(beta,
                              clust);
    }

    // initialize useful quantities
    k = beta.n_rows;
    probs.resize(k+1);
    probs.fill(0);

    // compute probabilities vector
    for(arma::uword j = 0; j < k; j++) {
      int nj = (int) arma::sum(clust == j);
      probs(j) = log(nj - sigma_PY) - log(2 * M_PI * sigma2) / 2 -
        (pow(y(i) - arma::dot(covs.row(i), beta.row(j)), 2) / (2 * sigma2));
    }

    prob_temp = 0.0;
    for(arma::uword l = 0; l < napprox; l++){
      prob_temp += arma::normpdf(y(i), arma::dot(covs.row(i), betatemp.row(l)), sqrt(sigma2));
    }
    probs(k) = log(mass + k * sigma_PY) + log(prob_temp / napprox);

    // sample new
    int temp_cl = rintnunif_log(probs);
    clust(i) = temp_cl;

    if(temp_cl == k){

      beta.resize(k+1, beta.n_cols);

      Sbn = arma::inv(arma::inv(Sb0) + arma::trans(covs.row(i)) * covs.row(i) / sigma2);
      betan = Sbn * (arma::inv(Sb0) * beta0 + (arma::trans(covs.row(i)) * y(i)) / sigma2);
      beta.row(k) = arma::trans(arma::mvnrnd(betan, Sbn));

    }
  }
}
