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

/*==================================================================================
 * update u - slice sampler
 *
 * args:
 * - clust:   vector of allocation
 * - w:       vector of stick weights
 * - u:       vector of uniform values
 *
 * Void function
 ==================================================================================*/

void update_u_SLI(arma::vec clust,
                  arma::vec w,
                  arma::vec &u){
  int nel = clust.n_elem;

  for(arma::uword el = 0; el < nel; el++){
    u[el] = arma::randu() * w(clust[el]);
  }
}

/*----------------------------------------------------------------------
 *
 * UNIVARIATE
 * LOCATION KERNEL
 * slice sampler functions
 *
 *----------------------------------------------------------------------
 */

/*==================================================================================
 * Accelerate - UNIVARIATE slice sampler - LOCATION
 * acceleration step for reshuffling the parameters
 * given an allocation
 *
 * args:
 * - data:  vector of observation
 * - mu:    vector of location component
 * - s2:    scale component
 * - v:     vector of stick break components
 * - w:     vector of stick weights
 * - clust: vector of allocation
 * - m0:    mean of location's prior distribution (scalar)
 * - s20:   prior variance of the location component
 * - a0:    shape of prior Gamma distribution on the scale component (scalar)
 * - b0:    rate of prior Gamma distribution on the scale component (scalar)
 * - mass:  mass of PY process
 * - sigma_PY: second parameter PY
 *
 * Void function
 ==================================================================================*/

void accelerate_SLI_PY_L(arma::vec data,
                         arma::vec &mu,
                         double &s2,
                         arma::vec &v,
                         arma::vec &w,
                         arma::vec clust,
                         double m0,
                         double s20,
                         double a0,
                         double b0,
                         double mass,
                         double sigma_PY){
  double s_temp, m_temp, xtemp, ytemp;
  double accu_temp = 0.0;

  for(arma::uword j = 0; j < mu.n_elem; j++){
    int nj  = arma::accu(clust == j);
    int ngj = arma::accu(clust > j);
    arma::vec tdata = data.elem(arma::find(clust == j));

    if(tdata.n_elem > 0){
      s_temp = 1 / (nj/s2 + 1/s20);
      m_temp = s_temp * (m0 / s20 + arma::accu(tdata)/s2);

      mu[j] = arma::randn() * sqrt(s_temp) + m_temp;
      accu_temp += arma::accu(pow(tdata - mu[j], 2));
    } else {
      s_temp = s20;
      m_temp = m0;
      mu[j] = arma::randn() * sqrt(s_temp) + m_temp;
    }

    xtemp = arma::randg(arma::distr_param(1.0 - sigma_PY + nj, 1.0));
    ytemp = arma::randg(arma::distr_param(mass + (j + 1) * sigma_PY + ngj, 1.0));
    v[j]  = xtemp / (xtemp + ytemp);

    if(j == 0){
      w[j] = v[j];
    }else{
      w[j] = v[j] * ((1 - v[j-1]) * w[j - 1]) / v[j-1];
    }
  }

  s2 = 1.0 / arma::randg(arma::distr_param((a0 + data.n_elem / 2), 1.0 / (b0 + accu_temp / 2)));
}

/*==================================================================================
 * Hyper-accelerate - UNIVARIATE slice sampler - LOCATION
 * acceleration step for reshuffling the parameters
 * given an allocation
 *
 * args:
 * - mu:    vector of location components
 * - clust: vector of allocations
 * - m0:    mean of location's prior distribution (scalar)
 * - s20:   prior variance of the location component
 * - m1:    hyperparameter, mean of distribution of mu
 * - k1:    hyperparameter, scale factor of distribution of mu
 * - a1:    hyperparameter, shape parameter of distribution of s20
 * - b1:    hyperparameter, scale parameter of distribution of s20
 *
 * Void function
 ==================================================================================*/

void hyper_accelerate_SLI_L(arma::vec mu,
                            arma::vec clust,
                            double &m0,
                            double &s20,
                            double m1,
                            double k1,
                            double a1,
                            double b1){

  double m_temp, k_temp, a_temp, b_temp, mu_m, acc = 0, acc2 = 0;
  int k = 0;

  for(arma::uword j = 0; j < mu.n_elem; j++){
    if(arma::accu(clust == j) != 0){
      acc  += mu(j);
      acc2 += pow(mu(j), 2);
      k += 1;
    }
  }

  mu_m = acc / k;
  k_temp = (k1 + k);
  m_temp = ((m1 * k1) + k * mu_m) / k_temp;
  a_temp = a1 + (k / 2.0);
  b_temp = b1 + ((acc2 - 2 * mu_m * acc + k * pow(mu_m, 2)) + (k * k1 * pow(mu_m - m1, 2)) / (k_temp)) / 2;

  s20 = 1.0 / arma::randg(arma::distr_param(a_temp, 1 / b_temp));
  m0  = arma::randn() * sqrt(s20 / k_temp) + m_temp;
}

/*==================================================================================
 * grow parameters - UNIVARIATE slice sampler - LOCATION
 * growing up the parameter vectors
 * till reaching the condition sum(w) > u_i, for all i
 *
 * args:
 * - mu:        vector, each element a mean
 * - v:         vector of stick break components
 * - w:         vector of stick weights
 * - u:         vector of uniform values
 * - m0:        mean of location's prior distribution (scalar)
 * - s20:       variance of the location component
 * - mass:      mass parameter
 * - n:         number of observations
 * - sigma_PY:  second parameter PY
 *
 * Void function
 ==================================================================================*/

void grow_param_SLI_PY_L(arma::vec &mu,
                         arma::vec &v,
                         arma::vec &w,
                         arma::vec u,
                         double m0,
                         double s20,
                         double mass,
                         int n,
                         double sigma_PY){
  double xtemp, ytemp;
  double w_sum = arma::accu(w);
  int new_val;
  int k_old = mu.n_elem;
  int k = w.n_elem;

  while(sum((1 - u) < w_sum) < n){

    k = w.n_elem;
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
    w_sum += w[k];
  }

  if(w.n_elem > k_old){
    new_val = w.n_elem - k_old;
    arma::vec mu_temp = arma::randn(new_val) * sqrt(s20) + m0;
    mu = arma::join_cols(mu, mu_temp);
  }
}

/*==================================================================================
 * grow parameters - UNIVARIATE slice sampler - LOCATION - INDEP
 * growing up the parameter vectors
 * till reaching the condition sum(w) > u_i, for all i
 *
 * args:
 * - mu:        vector, each element a mean
 * - v:         vector of stick break components
 * - w:         vector of stick weights
 * - xi:        weights of slice sampler
 * - u:         vector of uniform values
 * - m0:        mean of location's prior distribution (scalar)
 * - s20:       variance of the location component
 * - mass:      mass parameter
 * - n:         number of observations
 * - sigma_PY:  second parameter PY
 *
 * Void function
 ==================================================================================*/

void grow_param_indep_SLI_PY_L(arma::vec &mu,
                               arma::vec &v,
                               arma::vec &w,
                               arma::vec &xi,
                               arma::vec u,
                               double m0,
                               double s20,
                               double mass,
                               int n,
                               double sigma_PY,
                               double param_seq_one,
                               double param_seq_two){
  double xtemp, ytemp;
  double xi_sum = arma::accu(xi);
  int new_val;
  int k_old = mu.n_elem;
  int k = xi.n_elem;

  while(sum((1 - u) <= xi_sum) < n){

    k = xi.n_elem;
    v.resize(k + 1);
    w.resize(k + 1);
    xi.resize(k + 1);

    xtemp = arma::randg(arma::distr_param(1.0 - sigma_PY, 1.0));
    ytemp = arma::randg(arma::distr_param(mass + (k + 1) * sigma_PY, 1.0));
    v[k]  = xtemp / (xtemp + ytemp);

    if(k == 0){
      w[k] = v[k];
    }else{
      w[k] = v[k] * ((1 - v[k-1]) * w[k - 1]) / v[k-1];
    }
    xi[k] = xi[k - 1] * (param_seq_one + k * param_seq_two) / (param_seq_one + 1 + k * param_seq_two);
    xi_sum += xi[k];
  }

  if(xi.n_elem > k_old){
    new_val = xi.n_elem - k_old;
    arma::vec mu_temp = arma::randn(new_val) * sqrt(s20) + m0;
    mu = arma::join_cols(mu, mu_temp);
  }
}

/*==================================================================================
 * update cluster - UNIVARIATE slice sampler - LOCATION
 *
 * args:
 * - data:    vector of observation
 * - mu:      vector, each element a mean
 * - s2:      variance componet
 * - clust:   vector, each (integer) value is the cluster of the corresp. obs.
 * - w:       vector of stick weights
 * - u:       vector of uniform values
 *
 * Void function
 ==================================================================================*/

void update_cluster_SLI_L(arma::vec data,
                          arma::vec mu,
                          double s2,
                          arma::vec &clust,
                          arma::vec w,
                          arma::vec u){
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
    } else {
      probs.resize(index_use.n_elem);
      for(arma::uword j = 0; j < index_use.n_elem; j++){
        probs[j] = log(arma::normpdf(data[i], mu(index_use[j]), sqrt(s2)));
      }
      sampled = rintnunif_log(probs);
      clust[i] = index_use[sampled];
    }
  }
}

/*==================================================================================
 * update cluster - UNIVARIATE slice sampler - LOCATION - INDEP
 *
 * args:
 * - data:    vector of observation
 * - mu:      vector, each element a mean
 * - s2:      variance componet
 * - clust:   vector, each (integer) value is the cluster of the corresp. obs.
 * - w:       vector of stick weights
 * - xi:      weights for indep slice sampler
 * - u:       vector of uniform values
 *
 * Void function
 ==================================================================================*/

void update_cluster_indep_SLI_L(arma::vec data,
                                arma::vec mu,
                                double s2,
                                arma::vec &clust,
                                arma::vec w,
                                arma::vec xi,
                                arma::vec u){
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
      if(xi[r] > u[i]){
        siz++;
        index_use.resize(siz);
        index_use[siz - 1] = index[r];
      }
    }

    if(index_use.n_elem == 1){
      clust[i] = index_use[0];
    } else {
      probs.resize(index_use.n_elem);
      for(arma::uword j = 0; j < index_use.n_elem; j++){
        probs[j] = log(w(index_use[j])) - log(xi(index_use[j])) +
          log(arma::normpdf(data[i], mu(index_use[j]), sqrt(s2)));
      }
      sampled = rintnunif_log(probs);
      clust[i] = index_use[sampled];
    }
  }
}


/*----------------------------------------------------------------------
 *
 * UNIVARIATE
 * LOCATION-SCALE KERNEL
 * slice sampler functions
 *
 *----------------------------------------------------------------------
 */

/*==================================================================================
 * Accelerate - UNIVARIATE slice sampler - LOCATION SCALE
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
 * - m0:    mean of distribution of location components
 * - k0:    scale factor of distribution of location components
 * - a0:    shape of distribution on the scale component (scalar)
 * - b0:    scale of distribution on the scale component (scalar)
 * - mass:  mass of PY process
 * - sigma_PY: second parameter PY
 *
 * Void function
 ==================================================================================*/

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
  double xtemp, ytemp, kn, mn, an, bn, data_m;
  arma::vec tdata;

  for(arma::uword j = 0; j < mu.n_elem; j++){
    int nj  = arma::accu(clust == j);
    int ngj = arma::accu(clust > j);
    tdata = data.elem(arma::find(clust == j));
    data_m = sum(tdata) / nj;

    if(tdata.n_elem > 0){
      kn = (k0 + nj);
      mn = ((m0 * k0) + nj * data_m) / kn;
      an = a0 + (nj / 2.0);
      bn = b0 + (arma::accu(pow(tdata - data_m, 2)) + (nj * k0 * pow(data_m - m0, 2)) / (kn)) / 2;

      s2[j] = 1.0 / arma::randg(arma::distr_param(an, 1 / bn));
      mu[j] = arma::randn() * sqrt(s2[j] / kn) + mn;
    } else {
      s2[j] = 1.0 / arma::randg(arma::distr_param(a0, 1 / b0));
      mu[j] = arma::randn() * sqrt(s2[j] / k0) + m0;
    }

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

/*==================================================================================
 * Hyper-accelerate - UNIVARIATE slice sampler - LOCATION
 * acceleration step for reshuffling the parameters
 * given an allocation
 *
 * args:
 * - mu:    vector of location component
 * - s2:    vector of scale component
 * - m0:    mean of distribution of location components
 * - k0:    scale factor of distribution of location components
 * - a0:    shape of distribution on the scale component (scalar)
 * - b0:    scale of distribution on the scale component (scalar)
 * - m1:    hyperparameter, location component of m0
 * - s21:   hyperparameter, scale component of m0
 * - tau1:  hyperparameter, shape component of k0
 * - tau2:  hyperparameter, rate component of k0
 * - a1:    hyperparameter, shape component of b0
 * - b1:    hyperparameter, rate component of b0
 *
 * Void function
 ==================================================================================*/

void hyper_accelerate_SLI(arma::vec mu,
                          arma::vec s2,
                          arma::vec clust,
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

  double m_temp, s2_temp, a_temp, b_temp, tau1_temp, tau2_temp;
  double m2m0s2_acc = 0, ms2_acc = 0, s2_acc = 0;
  int k = 0;

  for(arma::uword j = 0; j < mu.n_elem; j++){
    if(arma::accu(clust == j) != 0){
      ms2_acc  += mu(j) / s2(j);
      m2m0s2_acc += pow(mu(j) - m0, 2) / s2(j);
      s2_acc   += 1 / s2(j);
      k += 1;
    }
  }

  tau1_temp = tau1 + k / 2;
  tau2_temp = tau2 + (m2m0s2_acc) / 2;
  k0 = arma::randg(arma::distr_param(tau1_temp, 1 / tau2_temp));

  s2_temp = 1 / ( 1 / s21 + k0 * s2_acc );
  m_temp  = s2_temp * ( m1 / s21 + k0 * ms2_acc );
  m0  = arma::randn() * sqrt(s2_temp) + m_temp;

  a_temp = a1 + k * a0;
  b_temp = b1 + s2_acc;
  b0 = arma::randg(arma::distr_param(a_temp, 1 / b_temp));
}

/*==================================================================================
 * grow parameters - UNIVARIATE slice sampler - PY - LOCATION SCALE
 * growing up the parameter vectors
 * till reaching the condition sum(w) > u_i, for all i
 *
 * args:
 * - mu:      vector, each element a mean
 * - s2:      vector, each element a variance
 * - v:       vector of stick break components
 * - w:       vector of stick weights
 * - u:       vector of uniform values
 * - m0:      mean of distribution of location components
 * - k0:      scale factor of distribution of location components
 * - a0:      shape of distribution on the scale component
 * - b0:      scale of distribution on the scale component
 * - mass:    mass parameter
 * - n:       number of observations
 * - sigma_PY: second parameter PY
 *
 * Void function
 ==================================================================================*/

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
                       double sigma_PY,
                       int &bound){
  double xtemp, ytemp;
  double w_sum = arma::accu(w);
  int new_val;
  int k_old = mu.n_elem;
  int k = w.n_elem;

  while(sum((1 - u) < w_sum) < n){

    if(k < pow(10.0, 5.0) - 1){
      k = w.n_elem;
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
    } else {
      bound += 1;
      k = w.n_elem;
      v.resize(k + 1);
      w.resize(k + 1);

      v[k]  = 1;

      if(k == 0){
        w[k] = v[k];
      }else{
        w[k] = v[k] * ((1 - v[k-1]) * w[k - 1]) / v[k-1];
      }
      w_sum = arma::accu(w);
    }

  }

  if(w.n_elem > k_old){
    new_val = w.n_elem - k_old;

    arma::vec s2_temp = 1.0 / arma::randg(new_val, arma::distr_param(a0, 1.0 / b0));
    arma::vec mu_temp = arma::randn(new_val) % sqrt(s2_temp / k0) + m0;

    mu = arma::join_cols(mu, mu_temp);
    s2 = arma::join_cols(s2, s2_temp);
  }
}

/*==================================================================================
 * grow parameters - UNIVARIATE slice sampler - PY - LOCATION SCALE - INDEP
 * growing up the parameter vectors
 * till reaching the condition sum(w) > u_i, for all i
 *
 * args:
 * - mu:      vector, each element a mean
 * - s2:      vector, each element a variance
 * - v:       vector of stick break components
 * - w:       vector of stick weights
 * - xi:      vector of weights for independent slice efficient
 * - u:       vector of uniform values
 * - m0:      mean of distribution of location components
 * - k0:      scale factor of distribution of location components
 * - a0:      shape of distribution on the scale component
 * - b0:      scale of distribution on the scale component
 * - mass:    mass parameter
 * - n:       number of observations
 * - sigma_PY: second parameter PY
 *
 * Void function
 ==================================================================================*/

void grow_param_indep_SLI_PY(arma::vec &mu,
                             arma::vec &s2,
                             arma::vec &v,
                             arma::vec &w,
                             arma::vec &xi,
                             arma::vec u,
                             double m0,
                             double k0,
                             double a0,
                             double b0,
                             double mass,
                             int n,
                             double sigma_PY,
                             double param_seq_one,
                             double param_seq_two,
                             int &bound){
  double xtemp, ytemp;
  double xi_sum = arma::accu(xi);
  int new_val;
  int k_old = mu.n_elem;
  int k = xi.n_elem;

  while(sum((1 - u) <= xi_sum) < n){

    if(k < pow(10.0, 5.0) - 1){
      k = xi.n_elem;
      v.resize(k + 1);
      w.resize(k + 1);
      xi.resize(k+1);

      xtemp = arma::randg(arma::distr_param(1.0 - sigma_PY, 1.0));
      ytemp = arma::randg(arma::distr_param(mass + (k + 1) * sigma_PY, 1.0));
      v[k]  = xtemp / (xtemp + ytemp);

      if(k == 0){
        w[k] = v[k];
      }else{
        w[k] = v[k] * ((1 - v[k-1]) * w[k - 1]) / v[k-1];
      }

      // xi[k] = xi[k - 1] * (param_seq_one + k * param_seq_two) / (param_seq_one + 1 + k * param_seq_two);
      // xi_sum += xi[k];

      xi[k] = xi[k - 1] * (mass + k * sigma_PY) / (mass + 1 + k * sigma_PY);
      xi_sum += xi[k];
    } else {

      bound += 1;
      double w_temp = arma::accu(w);
      k = xi.n_elem;
      v.resize(k + 1);
      w.resize(k + 1);
      xi.resize(k+1);

      v[k] = 1;
      w[k] = 1 - w_temp;

      xi[k] = 1 - xi_sum;
      xi_sum += xi[k];

      // xi[k] = xi[k - 1] * (mass + k * sigma_PY) / (mass + 1 + k * sigma_PY);
      // xi_sum += xi[k];
    }

    Rcpp::checkUserInterrupt();
    // Rcpp::Rcout << "\n\n" << k << "\n\n" << sum(log(1 - u) < log(xi_sum)) << "\n\n" << u.n_elem <<
    //   "\n\n" << sum(u <= 0) << "\n\n" << sum(u >= 1) << "\n\n"<< log(xi_sum) << "\n\n";
    //
    // if(k > 10000){
    //   Rcpp::Rcout << "\n\n" << u.t() << "\n\n" << (log(1 - u)).t() << "\n\n";
    // }
  }

  if(xi.n_elem > k_old){
    new_val = w.n_elem - k_old;

    arma::vec s2_temp = 1.0 / arma::randg(new_val, arma::distr_param(a0, 1.0 / b0));
    arma::vec mu_temp = arma::randn(new_val) % sqrt(s2_temp / k0) + m0;

    mu = arma::join_cols(mu, mu_temp);
    s2 = arma::join_cols(s2, s2_temp);
  }
}

/*==================================================================================
 * update cluster - UNIVARIATE conditional Slice sampler
 *
 * args:
 * - data:    vector of observation
 * - mu:      vector of location component
 * - s2:      vector of scale component
 * - clust:   vector of allocations
 * - w:       vector of stick weights
 * - u:       vector of uniform values
 *
 * Void function
 ==================================================================================*/

void update_cluster_SLI(arma::vec data,
                        arma::vec mu,
                        arma::vec s2,
                        arma::vec &clust,
                        arma::vec w,
                        arma::vec u){
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
    } else {
      probs.resize(index_use.n_elem);
      for(arma::uword j = 0; j < index_use.n_elem; j++){
        probs[j] = log(arma::normpdf(data[i], mu(index_use[j]), sqrt(s2(index_use[j]))));
      }
      sampled = rintnunif_log(probs);
      clust[i] = index_use[sampled];
    }
  }
}

/*==================================================================================
 * update cluster - UNIVARIATE conditional Slice sampler - INDEP
 *
 * args:
 * - data:    vector of observation
 * - mu:      vector of location component
 * - s2:      vector of scale component
 * - clust:   vector of allocations
 * - w:       vector of stick weights
 * - xi:      vector of weights for independent slice efficient
 * - u:       vector of uniform values
 *
 * Void function
 ==================================================================================*/

void update_cluster_indep_SLI(arma::vec data,
                              arma::vec mu,
                              arma::vec s2,
                              arma::vec &clust,
                              arma::vec w,
                              arma::vec xi,
                              arma::vec u){
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
      if(xi[r] > u[i]){
        siz++;
        index_use.resize(siz);
        index_use[siz - 1] = index[r];
      }
    }

    if(index_use.n_elem == 1){
      clust[i] = index_use[0];
    } else {
      probs.resize(index_use.n_elem);
      for(arma::uword j = 0; j < index_use.n_elem; j++){
        probs[j] = log(w(index_use[j])) - log(xi(index_use[j])) +
          log(arma::normpdf(data[i], mu(index_use[j]), sqrt(s2(index_use[j]))));
      }
      sampled = rintnunif_log(probs);
      clust[i] = index_use[sampled];
    }
  }
}


/*----------------------------------------------------------------------
 *
 * MULTIVARIATE
 * LOCATION KERNEL
 * slice sampler functions
 *
 *----------------------------------------------------------------------
 */

/*==================================================================================
 * Accelerate - MULTIVARIATE slice sampler - LOCATION
 * acceleration step for reshuffling the parameters
 * given an allocation
 *
 * args:
 * - data:  matrix of observations
 * - mu:    matrix of location components
 * - s2:    matrix - scale component
 * - v:     vector of stick break components
 * - w:     vector of stick weights
 * - clust: vector of allocation
 * - m0:    mean of location's prior distribution
 * - S20:   covariance matrix of the location component distribution
 * - S0:    matrix of Inverse Wishart distribution on s2
 * - n0:    degree of freedom of Inverse Wishart distribution on s2
 * - mass:  mass parameter
 * - sigma_PY: discount parameter
 *
 * Void function
 ==================================================================================*/

void accelerate_SLI_PY_mv_L(arma::mat data,
                            arma::mat &mu,
                            arma::mat &s2,
                            arma::vec &v,
                            arma::vec &w,
                            arma::vec clust,
                            arma::vec m0,
                            arma::mat S20,
                            arma::mat S0,
                            double n0,
                            double mass,
                            double sigma_PY){
  double xtemp, ytemp;
  arma::mat S2_temp = s2;
  arma::vec m_temp = m0;
  arma::mat tdata, cdata, S_update(s2);
  S_update.fill(0);
  int nj, ngj;
  int n = data.n_rows;

  // loop over the clusters
  for(arma::uword j = 0; j < mu.n_rows; j++){
    // initialize intra cluster quantities
    nj  = arma::accu(clust == j);
    ngj = arma::accu(clust > j);
    tdata = data.rows(arma::find(clust == j));

    // update parameters
    if(tdata.n_rows > 0){
      S2_temp = arma::inv(arma::inv(S20) + nj * arma::inv(s2));
      m_temp = S2_temp * (arma::inv(s2) * arma::trans(sum(tdata, 0)) + arma::inv(S20) * m0);
      mu.row(j)   = arma::trans(arma::mvnrnd(m_temp, S2_temp));

      cdata = tdata - arma::repmat(mu.row(j), nj, 1);
      S_update += arma::trans(cdata) * cdata;

    } else {
      S2_temp = S20;
      m_temp = m0;
      mu.row(j)   = arma::trans(arma::mvnrnd(m_temp, S2_temp));
    }

    // sample from the posterior distributions
    xtemp = arma::randg(arma::distr_param(1.0 - sigma_PY + nj, 1.0));
    ytemp = arma::randg(arma::distr_param(mass + (j + 1) * sigma_PY + ngj, 1.0));
    v[j]  = xtemp / (xtemp + ytemp);

    if(j == 0){
      w[j] = v[j];
    }else{
      w[j] = v[j] * ((1 - v[j-1]) * w[j - 1]) / v[j-1];
    }
  }

  arma::mat Sn = S0 + S_update;
  s2 = arma::inv(arma::wishrnd(arma::inv(Sn), (n0 + n)));
}

/*==================================================================================
 * Hyper-accelerate - MULTIVARIATE slice sampler - LOCATION
 *
 * args:
 * - mu:      vector of location component
 * - m0:      mean of location's prior distribution (scalar)
 * - clust:   vector of allocation
 * - S20:     variance of the location component
 * - m1:      hyperparameter, location component of m0
 * - k1:      hyperparameter, scale parameter for variance of m0
 * - theta1:  hyperparameter, df of S20 distribution
 * - Theta2:  hyperparameter, matrix of S20 distribution
 *
 * Void function
 ==================================================================================*/

void hyper_accelerate_SLI_mv_L(arma::mat mu,
                               arma::vec &m0,
                               arma::vec clust,
                               arma::mat &S20,
                               arma::vec m1,
                               double k1,
                               double theta1,
                               arma::mat Theta1){
  arma::mat tmu(mu);
  int k = 0;
  for(arma::uword j = 0; j < mu.n_rows; j++){
    if(arma::accu(clust == j) != 0){
      tmu.row(k) = mu.row(j);
      k++;
    }
  }
  tmu.resize(k, tmu.n_cols);

  // compute the updated parameters
  arma::vec mum = arma::trans(mean(tmu, 0));
  double k1n = k1 + k;
  arma::vec m1n = ((m1 * k1) + arma::trans(sum(tmu, 0)))/k1n;
  arma::mat cmu = tmu - arma::repmat(mum.t(), k, 1);
  double theta1n = theta1 + k;
  arma::mat Theta1n = Theta1 + arma::trans(cmu) * cmu + ((k1*k)/k1n) * (mum - m0) * arma::trans(mum - m0);

  // sample from the posterior distributions
  S20 = arma::inv(arma::wishrnd(arma::inv(Theta1n), theta1n));
  m0  = arma::mvnrnd(m1n, S20/k1n);
}

/*==================================================================================
 * grow parameters - MULTIVARIATE slice sampler - LOCATION
 * growing up the parameter vectors
 * till reaching the condition sum(w) > u_i, for all i
 *
 * args:
 * - mu:      matrix, each row a mean
 * - s2:      matrix, covariance
 * - v:       vector of stick break components
 * - w:       vector of stick weights
 * - u:       vector of uniform values
 * - m0:      mean of location's prior distribution
 * - S20:     covariance matrix of the location component distribution
 * - mass:    mass parameter
 * - n:       number of observations
 * - sigma_PY: discount parameter
 * - max_len: maximum span
 * - n_over:  number of times reached max_len
 *
 * Void function
 ==================================================================================*/

void grow_param_SLI_PY_mv_L(arma::mat &mu,
                            arma::vec &v,
                            arma::vec &w,
                            arma::vec u,
                            arma::vec m0,
                            arma::mat S20,
                            double mass,
                            int n,
                            double sigma_PY){

  double xtemp, ytemp;
  double w_sum = arma::accu(w);
  int k_old = mu.n_rows;
  int k = w.n_elem;
  int k_new;

  while(sum((1 - u) < w_sum) < n){

    k = w.n_elem;
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

  k_new = w.n_elem;
  mu.resize(k_new, mu.n_cols);
  for(arma::uword j = k_old; j < k_new; j++){
    mu.row(j) = arma::trans(arma::mvnrnd(m0, S20));
  }
}

/*==================================================================================
 * grow parameters - MULTIVARIATE slice sampler - LOCATION - INDEP
 * growing up the parameter vectors
 * till reaching the condition sum(w) > u_i, for all i
 *
 * args:
 * - mu:      matrix, each row a mean
 * - s2:      matrix, covariance
 * - v:       vector of stick break components
 * - w:       vector of stick weights
 * - u:       vector of uniform values
 * - m0:      mean of location's prior distribution
 * - S20:     covariance matrix of the location component distribution
 * - mass:    mass parameter
 * - n:       number of observations
 * - sigma_PY: discount parameter
 * - max_len: maximum span
 * - n_over:  number of times reached max_len
 *
 * Void function
 ==================================================================================*/

void grow_param_indep_SLI_PY_mv_L(arma::mat &mu,
                                  arma::vec &v,
                                  arma::vec &w,
                                  arma::vec &xi,
                                  arma::vec u,
                                  arma::vec m0,
                                  arma::mat S20,
                                  double mass,
                                  int n,
                                  double sigma_PY,
                                  double param_seq_one,
                                  double param_seq_two){

  double xtemp, ytemp;
  double xi_sum = arma::accu(xi);
  int k_old = mu.n_rows;
  int k = w.n_elem;
  int k_new;

  while(sum((1 - u) <= xi_sum) < n){

    k = w.n_elem;
    v.resize(k+1);
    w.resize(k+1);
    xi.resize(k+1);

    xtemp = arma::randg(arma::distr_param(1.0 - sigma_PY, 1.0));
    ytemp = arma::randg(arma::distr_param(mass + (k + 1) * sigma_PY, 1.0));
    v[k]  = xtemp / (xtemp + ytemp);

    if(k == 0){
      w[k] = v[k];
    }else{
      w[k] = v[k] * ((1 - v[k-1]) * w[k - 1]) / v[k-1];
    }

    xi[k] = xi[k - 1] * (param_seq_one + k * param_seq_two) / (param_seq_one + 1 + k * param_seq_two);
    xi_sum += xi[k];
  }

  k_new = w.n_elem;
  mu.resize(k_new, mu.n_cols);
  for(arma::uword j = k_old; j < k_new; j++){
    mu.row(j) = arma::trans(arma::mvnrnd(m0, S20));
  }
}

/*==================================================================================
 * update cluster - MULTIVARIATE slice sampler - LOCATION
 *
 * args:
 * - data:    matrix of observation
 * - mu:      matrix of location component
 * - s2:      covariance matrix
 * - clust:   vector of allocations
 * - w:       vector of stick weights
 * - u:       vector of uniform values
 *
 * Void function
 ==================================================================================*/

void update_cluster_SLI_mv_L(arma::mat data,
                             arma::mat mu,
                             arma::mat s2,
                             arma::vec &clust,
                             arma::vec w,
                             arma::vec u){
  int n = data.n_rows;
  int k = mu.n_rows;
  int d = data.n_cols;
  arma::uvec index_use;
  arma::uvec index = arma::regspace<arma::uvec>(0, k - 1);
  arma::vec probs;
  int siz;
  int sampled;

  arma::mat rooti  = arma::trans(arma::inv(trimatu(arma::chol(s2))));
  double rsum = arma::sum(log(rooti.diag()));

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
    } else {
      probs.resize(index_use.n_elem);
      for(arma::uword j = 0; j < index_use.n_elem; j++){
        arma::vec cdata  = rooti * arma::trans(data.row(i) - mu.row(index_use[j])) ;
        probs[j]         = - (d / 2.0) * log(2.0 * M_PI) - 0.5 * arma::sum(cdata%cdata) + rsum;
      }
      sampled = rintnunif_log(probs);
      clust[i] = index_use[sampled];
    }
  }
}

/*==================================================================================
 * update cluster - MULTIVARIATE slice sampler - LOCATION - INDEP
 *
 * args:
 * - data:    matrix of observation
 * - mu:      matrix of location component
 * - s2:      covariance matrix
 * - clust:   vector of allocations
 * - w:       vector of stick weights
 * - u:       vector of uniform values
 *
 * Void function
 ==================================================================================*/

void update_cluster_indep_SLI_mv_L(arma::mat data,
                                   arma::mat mu,
                                   arma::mat s2,
                                   arma::vec &clust,
                                   arma::vec w,
                                   arma::vec xi,
                                   arma::vec u){
  int n = data.n_rows;
  int k = mu.n_rows;
  int d = data.n_cols;
  arma::uvec index_use;
  arma::uvec index = arma::regspace<arma::uvec>(0, k - 1);
  arma::vec probs;
  int siz;
  int sampled;

  arma::mat rooti  = arma::trans(arma::inv(trimatu(arma::chol(s2))));
  double rsum = arma::sum(log(rooti.diag()));

  for(arma::uword i = 0; i < n; i++){
    siz = 0;
    index_use.resize(1);
    for(arma::uword r = 0; r < k; r++){
      if(xi[r] > u[i]){
        siz++;
        index_use.resize(siz);
        index_use[siz - 1] = index[r];
      }
    }

    if(index_use.n_elem == 1){
      clust[i] = index_use[0];
    } else {
      probs.resize(index_use.n_elem);
      for(arma::uword j = 0; j < index_use.n_elem; j++){
        arma::vec cdata  = rooti * arma::trans(data.row(i) - mu.row(index_use[j])) ;
        probs[j]         = log(w(index_use[j])) - log(xi(index_use[j])) - (d / 2.0) * log(2.0 * M_PI) - 0.5 * arma::sum(cdata%cdata) + rsum;
      }
      sampled = rintnunif_log(probs);
      clust[i] = index_use[sampled];
    }
  }
}

/*----------------------------------------------------------------------
 *
 * MULTIVARIATE
 * LOCATION-SCALE KERNEL
 * slice sampler functions
 *
 *----------------------------------------------------------------------
 */

/*==================================================================================
 * Accelerate - MULTIVARIATE conditional Slice sampler - PY - LOCATION SCALE
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
 * - mass:  mass parameter
 * - sigma_PY: discount parameter
 *
 * Void function
 ==================================================================================*/

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
  arma::vec mn;
  arma::mat cdata, tdata, Sn;
  double nn, kn, xtemp, ytemp;
  int nj, ngj;

  // loop over the clusters
  for(arma::uword j = 0; j < mu.n_rows; j++){
    // initialize intra cluster quantities
    nj  = arma::accu(clust == j);
    ngj = arma::accu(clust > j);
    tdata = data.rows(arma::find(clust == j));

    // update parameters
    if(tdata.n_rows > 0){
      kn = k0 + nj;
      mn = ((m0 * k0) + arma::trans(sum(tdata, 0)))/kn;
      cdata = tdata - arma::repmat(mean(tdata, 0), nj, 1);
      nn = n0 + nj;
      Sn = S0 + arma::trans(cdata) * cdata + ((k0*nj)/kn) *
        (arma::trans(mean(tdata, 0)) - m0) * arma::trans(arma::trans(mean(tdata, 0)) - m0);
    } else {
      kn = k0;
      mn = m0;
      nn = n0;
      Sn = S0;
    }

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

/*==================================================================================
 * Hyper-accelerate - MULTIVARIATE slice sampler - LOCATION-SCALE
 *
 * args:
 * - mu:      vector of location component
 * - s2:      covariance matricies
 * - clust:   vector of allocations
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

void hyper_accelerate_SLI_mv_LS(arma::mat mu,
                                arma::cube s2,
                                arma::vec clust,
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
  int k = 0;

  for(arma::uword j = 0; j < mu.n_rows; j++){
    if(arma::accu(clust == j) != 0){
      accu_Sinv   += arma::inv(s2.slice(j));
      accu_Smu    += arma::inv(s2.slice(j)) * mu.row(j).t();
      accu_qform  += arma::as_scalar(mu.row(j) * arma::inv(s2.slice(j)) * mu.row(j).t());
      k ++;
    }
  }

  // update m0
  arma::mat S1n = arma::inv(arma::inv(S1) + k0 * accu_Sinv);
  arma::vec m1n = S1n * (arma::inv(S1) * m1 + k0 * accu_Smu);
  m0 = arma::mvnrnd(m1n, S1n);

  // update S0
  double theta1n = theta1 + k * n0;
  arma::mat Theta1n = arma::inv(arma::inv(Theta1) + accu_Sinv);
  S0 = arma::wishrnd(Theta1n, theta1n);

  // update k0
  double tau1n = tau1 + k;
  double tau2n = tau2 + accu_qform;
  k0 = arma::randg(arma::distr_param(tau1n, 1 / tau2n));
}

/*==================================================================================
 * grow parameters - MULTIVARIATE slice sampler - LOCATION SCALE
 * growing up the parameter vectors
 * till reaching the condition sum(w) > u_i, for all i
 *
 * args:
 * - mu:        matrix, each row a mean
 * - s2:        cube, each slice a covariance matrix
 * - v:         vector of stick break components
 * - w:         vector of stick weights
 * - u:         vector of uniform values
 * - m0:        mean of location distribution
 * - k0:        scale factor of location distribution
 * - S20:       matrix of the scale distribution
 * - n0:        df of scale distribution
 * - mass:      mass parameter
 * - n:         number of observations
 * - sigma_PY:  discount parameter
 * - max_len:   maximum span
 * - n_over:    number of times reached max_len
 *
 * Void function
 ==================================================================================*/

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
  double xtemp, ytemp;
  double w_sum = arma::accu(w);
  int k_old = mu.n_rows;
  int k = w.n_elem;
  int k_new;

  while(sum((1 - u) < w_sum) < n){

    k = w.n_elem;
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

  k_new = w.n_elem;
  mu.resize(k_new, mu.n_cols);
  s2.resize(s2.n_rows, s2.n_cols, k_new);

  for(arma::uword j = k_old; j < k_new; j++){
    s2.slice(j) = arma::inv(arma::wishrnd(arma::inv(S0), n0));
    mu.row(j) = arma::trans(arma::mvnrnd(m0, s2.slice(j)/k0));
  }
}

/*==================================================================================
 * grow parameters - MULTIVARIATE slice sampler - LOCATION SCALE - INDEP
 * growing up the parameter vectors
 * till reaching the condition sum(w) > u_i, for all i
 *
 * args:
 * - mu:        matrix, each row a mean
 * - s2:        cube, each slice a covariance matrix
 * - v:         vector of stick break components
 * - w:         vector of stick weights
 * - u:         vector of uniform values
 * - m0:        mean of location distribution
 * - k0:        scale factor of location distribution
 * - S20:       matrix of the scale distribution
 * - n0:        df of scale distribution
 * - mass:      mass parameter
 * - n:         number of observations
 * - sigma_PY:  discount parameter
 * - max_len:   maximum span
 * - n_over:    number of times reached max_len
 *
 * Void function
 ==================================================================================*/

void grow_param_indep_SLI_PY_mv(arma::mat &mu,
                                arma::cube &s2,
                                arma::vec &v,
                                arma::vec &w,
                                arma::vec &xi,
                                arma::vec u,
                                arma::vec m0,
                                double k0,
                                arma::mat S0,
                                double n0,
                                double mass,
                                int n,
                                double sigma_PY,
                                double param_seq_one,
                                double param_seq_two){
  double xtemp, ytemp;
  double xi_sum = arma::accu(xi);
  int k_old = mu.n_rows;
  int k = w.n_elem;
  int k_new;

  while(sum((1 - u) <= xi_sum) < n){

    k = w.n_elem;
    v.resize(k+1);
    w.resize(k+1);
    xi.resize(k+1);

    xtemp = arma::randg(arma::distr_param(1.0 - sigma_PY, 1.0));
    ytemp = arma::randg(arma::distr_param(mass + (k + 1) * sigma_PY, 1.0));
    v[k]  = xtemp / (xtemp + ytemp);

    if(k == 0){
      w[k] = v[k];
    }else{
      w[k] = v[k] * ((1 - v[k-1]) * w[k - 1]) / v[k-1];
    }

    xi[k] = xi[k - 1] * (param_seq_one + k * param_seq_two) / (param_seq_one + 1 + k * param_seq_two);
    xi_sum += xi[k];
  }

  k_new = w.n_elem;
  mu.resize(k_new, mu.n_cols);
  s2.resize(s2.n_rows, s2.n_cols, k_new);

  for(arma::uword j = k_old; j < k_new; j++){
    s2.slice(j) = arma::inv(arma::wishrnd(arma::inv(S0), n0));
    mu.row(j) = arma::trans(arma::mvnrnd(m0, s2.slice(j)/k0));
  }
}

/*==================================================================================
 * update cluster - MULTIVARIATE conditional Slice sampler
 *
 * args:
 * - data:    matrix of observation
 * - mu:      matrix of location component
 * - s2:      cube of scale component
 * - clust:   vector of allocations
 * - w:       vector of stick weights
 * - u:       vector of uniform values
 *
 * Void function
 ==================================================================================*/

void update_cluster_SLI_mv(arma::mat data,
                           arma::mat mu,
                           arma::cube s2,
                           arma::vec &clust,
                           arma::vec w,
                           arma::vec u){
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
    } else {
      probs.resize(index_use.n_elem);
      for(arma::uword j = 0; j < index_use.n_elem; j++){
        arma::mat rooti  = arma::trans(arma::inv(trimatu(arma::chol(s2.slice(index_use[j])))));
        arma::vec cdata  = rooti * arma::trans(data.row(i) - mu.row(index_use[j])) ;
        probs[j]         = - (d / 2.0) * log(2.0 * M_PI) - 0.5 * arma::sum(cdata%cdata) + arma::sum(log(rooti.diag()));
      }
      sampled = rintnunif_log(probs);
      clust[i] = index_use[sampled];
    }
  }
}

/*==================================================================================
 * update cluster - MULTIVARIATE conditional Slice sampler - INDEP
 *
 * args:
 * - data:    matrix of observation
 * - mu:      matrix of location component
 * - s2:      cube of scale component
 * - clust:   vector of allocations
 * - w:       vector of stick weights
 * - u:       vector of uniform values
 *
 * Void function
 ==================================================================================*/

void update_cluster_indep_SLI_mv(arma::mat data,
                                 arma::mat mu,
                                 arma::cube s2,
                                 arma::vec &clust,
                                 arma::vec w,
                                 arma::vec xi,
                                 arma::vec u){
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
      if(xi[r] > u[i]){
        siz++;
        index_use.resize(siz);
        index_use[siz - 1] = index[r];
      }
    }

    if(index_use.n_elem == 1){
      clust[i] = index_use[0];
    } else {
      probs.resize(index_use.n_elem);
      for(arma::uword j = 0; j < index_use.n_elem; j++){
        arma::mat rooti  = arma::trans(arma::inv(trimatu(arma::chol(s2.slice(index_use[j])))));
        arma::vec cdata  = rooti * arma::trans(data.row(i) - mu.row(index_use[j])) ;
        probs[j]         = log(w(index_use[j])) - log(xi(index_use[j])) - (d / 2.0) * log(2.0 * M_PI) - 0.5 * arma::sum(cdata%cdata) + arma::sum(log(rooti.diag()));
      }
      sampled = rintnunif_log(probs);
      clust[i] = index_use[sampled];
    }
  }
}

/*----------------------------------------------------------------------
 *
 * MULTIVARIATE
 * PRODUCT KERNEL
 * slice sampler functions
 *
 *----------------------------------------------------------------------
 */

/*==================================================================================
 * Accelerate - MULTIVARIATE slice sampler - LOCATION SCALE
 * acceleration step for reshuffling the parameters
 * given an allocation
 *
 * args:
 * - data:      matrix of observation
 * - mu:        matrix of location component
 * - s2:        cube of scale component
 * - v:         vector of stick break components
 * - w:         vector of stick weights
 * - clust:     vector of allocation
 * - m0:        vector of location's prior distribution
 * - k0:        vector of location's variance tuning parameter, one for each dimension
 * - a0:        vector of shape parameters for the scale component, one for each dimension
 * - b0:        vector of scale parameters for the scale component, one for each dimension
 * - mass:      mass parameter
 * - sigma_PY:  discount parameter
 *
 * Void function
 ==================================================================================*/

void accelerate_SLI_PY_mv_P(arma::mat data,
                            arma::mat &mu,
                            arma::mat &s2,
                            arma::vec &v,
                            arma::vec &w,
                            arma::vec clust,
                            arma::vec m0,
                            arma::vec k0,
                            arma::vec a0,
                            arma::vec b0,
                            double mass,
                            double sigma_PY){
  arma::mat cdata, tdata;
  double kn, xtemp, ytemp, an, bn, mn, data_m;
  int nj, ngj;

  // loop over the clusters
  for(arma::uword j = 0; j < mu.n_rows; j++){

    // INITIALIZE CLUSTER RELATED QUANTITIES
    nj  = arma::accu(clust == j);
    ngj = arma::accu(clust > j);
    tdata = data.rows(arma::find(clust == j));

    // UPDATE THE LOCATIONS MARGINALLY
    if(tdata.n_rows > 0){
      for(arma::uword l = 0; l < mu.n_cols; l++){

        data_m = arma::accu(tdata.col(l)) / nj;
        kn = k0(l) + nj;
        mn = ((m0(l) * k0(l)) + nj * data_m) / kn;
        an = a0(l) + (nj / 2.0);
        bn = b0(l) + (arma::accu(pow(tdata.col(l) - data_m, 2)) + (nj * k0(l) * pow(data_m - m0(l), 2)) / (kn)) / 2;

        s2(j,l) = 1.0 / arma::randg(arma::distr_param(an, 1.0 / bn));
        mu(j,l) = arma::randn() * sqrt(s2(j,l) / kn) + mn;
      }
    } else {
      for(arma::uword k = 0; k < mu.n_cols; k++){
        s2(j,k) = 1.0 / arma::randg(arma::distr_param(a0(k), 1.0 / b0(k)));
        mu(j,k) = arma::randn() * sqrt(s2(j,k) / k0(k)) + m0(k);
      }
    }

    // SAMPLE THE WEIGHTS OF STICK-BREAKING REPRESENTATION
    xtemp = arma::randg(arma::distr_param(1.0 - sigma_PY + nj, 1.0));
    ytemp = arma::randg(arma::distr_param(mass + (j + 1) * sigma_PY + ngj, 1.0));
    v[j]  = xtemp / (xtemp + ytemp);

    if(j != 0){
      w[j] = v[j] * ((1 - v[j-1]) * w[j - 1]) / v[j-1];
    }else{
      w[j] = v[j];
    }
  }
}

/*==================================================================================
 * Hyper-accelerate - MULTIVARIATE slice sampler - PRODUCT
 * acceleration step for reshuffling the parameters
 * given an allocation
 *
 * args:
 * - mu:    vector of location component
 * - s2:    vector of scale component
 * - clust: vector of allocations
 * - m0:    mean of location's prior distribution (scalar)
 * - k0:    tuning parameter of variance of the location component
 * - a0:    shape parameters of scale component
 * - b0:    scale parameters of scale component
 * - m1:    hyperparameter, location component of m0
 * - s21:   hyperparameter, scale component of m0
 * - tau1:  hyperparameter, shape component of k0
 * - tau2:  hyperparameter, rate component of k0
 * - a1:    hyperparameter, shape component of b0
 * - b1:    hyperparameter, rate component of b0
 *
 * Void function
 ==================================================================================*/

void hyper_accelerate_SLI_mv_P(arma::mat mu,
                               arma::mat s2,
                               arma::vec clust,
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

  double m_temp, s2_temp, tau1_temp, tau2_temp, a_temp, b_temp;
  int k = 0;
  arma::vec ms2_acc(m0.n_elem, arma::fill::zeros);
  arma::vec m2m0s2_acc(m0.n_elem, arma::fill::zeros);
  arma::vec s2_acc(m0.n_elem, arma::fill::zeros);

  for(arma::uword l = 0; l < mu.n_elem; l++){
    if(arma::accu(clust == l) != 0){
      ms2_acc    += arma::trans(mu.row(l) / s2.row(l));
      m2m0s2_acc += arma::trans(pow(mu.row(l) - m0.t(), 2) / s2.row(l));
      s2_acc     += arma::trans(1 / s2.row(l));
      k += 1;
    }
  }

  for(arma::uword j = 0; j < mu.n_cols; j++){

    tau1_temp = tau1(j) + k / 2;
    tau2_temp = tau2(j) + m2m0s2_acc(j) / 2;
    k0(j) = arma::randg(arma::distr_param(tau1_temp, 1 / tau2_temp));

    s2_temp = 1 / ( 1 / s21(j) + k0(j) * s2_acc(j) );
    m_temp  = s2_temp * ( m1(j) / s21(j) + k0(j) * ms2_acc(j) );
    m0(j)  = arma::randn() * sqrt(s2_temp) + m_temp;

    a_temp = a1(j) + k * a0(j);
    b_temp = b1(j) + s2_acc(j);
    b0(j) = arma::randg(arma::distr_param(a_temp, 1 / b_temp));
  }
}

/*==================================================================================
 * grow parameters - MULTIVARIATE slice sampler - LOCATION SCALE
 * growing up the parameter vectors
 * till reaching the condition sum(w) > u_i, for all i
 *
 * args:
 * - mu:        matrix, each row a mean
 * - s2:        cube, each slice a covariance matrix
 * - v:         vector of stick break components
 * - w:         vector of stick weights
 * - u:         vector of uniform values
 * - m0:        vector of location's prior distribution
 * - k0:        vector of location's variance tuning parameter, one for each dimension
 * - a0:        vector of parameters for the scale component, one for each dimension
 * - b0:        vector of parameters for the scale component, one for each dimension
 * - mass:      mass parameter
 * - n:         number of observations
 * - sigma_PY:  discount parameter
 * - max_len:   maximum span
 * - n_over:    number of times reached max_len
 *
 * Void function
 ==================================================================================*/

void grow_param_SLI_PY_mv_P(arma::mat &mu,
                            arma::mat &s2,
                            arma::vec &v,
                            arma::vec &w,
                            arma::vec u,
                            arma::vec m0,
                            arma::vec k0,
                            arma::vec a0,
                            arma::vec b0,
                            double mass,
                            int n,
                            double sigma_PY){

  double xtemp, ytemp;
  double w_sum = arma::accu(w);
  int k_old = mu.n_rows;
  int k = w.n_elem;
  int k_new;

  while(sum((1 - u) < w_sum) < n){

    k = w.n_elem;
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

  k_new = w.n_elem;
  mu.resize(k_new, mu.n_cols);
  s2.resize(k_new, s2.n_cols);

  for(arma::uword j = k_old; j < k_new; j++){
    for(arma::uword k = 0; k < mu.n_cols; k++){
      s2(j,k) = 1.0 / arma::randg(arma::distr_param(a0(k), 1.0 / b0(k)));
      mu(j,k) = arma::randn() * sqrt(k0(k) * s2(j,k)) + m0(k);
    }
  }
}

/*==================================================================================
 * grow parameters - MULTIVARIATE slice sampler - LOCATION SCALE - INDEP
 * growing up the parameter vectors
 * till reaching the condition sum(w) > u_i, for all i
 *
 * args:
 * - mu:        matrix, each row a mean
 * - s2:        cube, each slice a covariance matrix
 * - v:         vector of stick break components
 * - w:         vector of stick weights
 * - u:         vector of uniform values
 * - m0:        vector of location's prior distribution
 * - k0:        vector of location's variance tuning parameter, one for each dimension
 * - a0:        vector of parameters for the scale component, one for each dimension
 * - b0:        vector of parameters for the scale component, one for each dimension
 * - mass:      mass parameter
 * - n:         number of observations
 * - sigma_PY:  discount parameter
 * - max_len:   maximum span
 * - n_over:    number of times reached max_len
 *
 * Void function
 ==================================================================================*/

void grow_param_indep_SLI_PY_mv_P(arma::mat &mu,
                                  arma::mat &s2,
                                  arma::vec &v,
                                  arma::vec &w,
                                  arma::vec &xi,
                                  arma::vec u,
                                  arma::vec m0,
                                  arma::vec k0,
                                  arma::vec a0,
                                  arma::vec b0,
                                  double mass,
                                  int n,
                                  double sigma_PY,
                                  double param_seq_one,
                                  double param_seq_two){

  double xtemp, ytemp;
  double xi_sum = arma::accu(xi);
  int k_old = mu.n_rows;
  int k = w.n_elem;
  int k_new;

  while(sum((1 - u) <= xi_sum) < n){

    k = w.n_elem;
    v.resize(k+1);
    w.resize(k+1);
    xi.resize(k+1);

    xtemp = arma::randg(arma::distr_param(1.0 - sigma_PY, 1.0));
    ytemp = arma::randg(arma::distr_param(mass + (k + 1) * sigma_PY, 1.0));
    v[k]  = xtemp / (xtemp + ytemp);

    if(k == 0){
      w[k] = v[k];
    }else{
      w[k] = v[k] * ((1 - v[k-1]) * w[k - 1]) / v[k-1];
    }

    xi[k] = xi[k - 1] * (param_seq_one + k * param_seq_two) / (param_seq_one + 1 + k * param_seq_two);
    xi_sum += xi[k];
  }

  k_new = w.n_elem;
  mu.resize(k_new, mu.n_cols);
  s2.resize(k_new, s2.n_cols);

  for(arma::uword j = k_old; j < k_new; j++){
    for(arma::uword k = 0; k < mu.n_cols; k++){
      s2(j,k) = 1.0 / arma::randg(arma::distr_param(a0(k), 1.0 / b0(k)));
      mu(j,k) = arma::randn() * sqrt(k0(k) * s2(j,k)) + m0(k);
    }
  }
}

/*==================================================================================
 * update cluster - MULTIVARIATE conditional Slice sampler
 *
 * args:
 * - data:    matrix of observation
 * - mu:      matrix of location component
 * - s2:      matrix of scale component
 * - clust:   vector, each (integer) value is the cluster of the corresp. obs.
 * - w:       vector of stick weights
 * - u:       vector of uniform values
 * - max_val: vector, number of already existent atoms
 * - iter:    current iteration
 * - new_val: vector of new values
 *
 * Void function
 ==================================================================================*/

void update_cluster_SLI_mv_P(arma::mat data,
                             arma::mat mu,
                             arma::mat s2,
                             arma::vec &clust,
                             arma::vec w,
                             arma::vec u){
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
    } else {
      probs.resize(index_use.n_elem);
      for(arma::uword j = 0; j < index_use.n_elem; j++){
        probs(j) = 0;
        for(arma::uword l = 0; l < d; l++){
          probs(j) +=  - log(2 * M_PI * s2(index_use(j),l)) / 2 - (pow(data(i,l) - mu(index_use(j),l), 2) / (2 * s2(index_use(j),l)));
        }
      }

      sampled = rintnunif(exp(probs));
      clust[i] = index_use[sampled];
    }
  }
}

/*==================================================================================
 * update cluster - MULTIVARIATE conditional Slice sampler - INDEP
 *
 * args:
 * - data:    matrix of observation
 * - mu:      matrix of location component
 * - s2:      matrix of scale component
 * - clust:   vector, each (integer) value is the cluster of the corresp. obs.
 * - w:       vector of stick weights
 * - u:       vector of uniform values
 * - max_val: vector, number of already existent atoms
 * - iter:    current iteration
 * - new_val: vector of new values
 *
 * Void function
 ==================================================================================*/

void update_cluster_indep_SLI_mv_P(arma::mat data,
                                   arma::mat mu,
                                   arma::mat s2,
                                   arma::vec &clust,
                                   arma::vec w,
                                   arma::vec xi,
                                   arma::vec u){
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
      if(xi[r] > u[i]){
        siz++;
        index_use.resize(siz);
        index_use[siz - 1] = index[r];
      }
    }

    if(index_use.n_elem == 1){
      clust[i] = index_use[0];
    } else {
      probs.resize(index_use.n_elem);
      for(arma::uword j = 0; j < index_use.n_elem; j++){
        probs(j) = 0;
        for(arma::uword l = 0; l < d; l++){
          probs(j) += log(w(index_use[j])) - log(xi(index_use[j])) - log(2 * M_PI * s2(index_use(j),l)) / 2 - (pow(data(i,l) - mu(index_use(j),l), 2) / (2 * s2(index_use(j),l)));
        }
      }

      sampled = rintnunif(exp(probs));
      clust[i] = index_use[sampled];
    }
  }
}


/*----------------------------------------------------------------------
 *
 * MIXTURE OF REGRESSION KERNELS
 * LOCATION-SCALE KERNEL
 * SLI functions
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

void accelerate_SLI_mv_MRK(arma::vec y,
                           arma::mat covs,
                           arma::mat &beta,
                           arma::vec &sigma2,
                           arma::vec &v,
                           arma::vec &w,
                           arma::vec clust,
                           arma::vec beta0,
                           arma::mat Sb0,
                           double a0,
                           double b0,
                           double mass,
                           double sigma_PY){

  arma::mat cdata;
  arma::mat tdata;
  arma::vec ty;

  double an, bn, xtemp, ytemp;
  arma::mat tSb;
  arma::vec tbeta0;

  int nj, ngj;

  // loop over the clusters
  for(arma::uword j = 0; j < beta.n_rows; j++){

    // INITIALIZE CLUSTER RELATED QUANTITIES
    nj  = arma::accu(clust == j);
    ngj = arma::accu(clust > j);
    tdata = covs.rows(arma::find(clust == j));
    ty = y.elem(arma::find(clust == j));

    // UPDATE THE LOCATIONS MARGINALLY
    if(tdata.n_rows > 0){
      cdata = ty - (tdata * beta.row(j).t());

      // update the variance
      an = a0 + nj/2;
      bn = b0 + (arma::accu(pow(cdata, 2))) / 2;
      sigma2(j) = 1.0 / arma::randg(arma::distr_param(an, 1.0 / (bn)));

      // update the coefficients
      tSb = arma::inv(arma::inv(Sb0) + arma::trans(tdata) * tdata / sigma2(j));
      tbeta0 = tSb * (arma::inv(Sb0) * beta0 + (arma::trans(tdata) * ty) / sigma2(j));;
      beta.row(j) = arma::trans(arma::mvnrnd(tbeta0, tSb));

    } else {

      sigma2(j) = 1.0 / arma::randg(arma::distr_param(a0, 1.0 / (b0)));
      beta.row(j) = arma::trans(arma::mvnrnd(beta0, Sb0));
    }

    // SAMPLE THE WEIGHTS OF STICK-BREAKING REPRESENTATION
    xtemp = arma::randg(arma::distr_param(1.0 - sigma_PY + nj, 1.0));
    ytemp = arma::randg(arma::distr_param(mass + (j + 1) * sigma_PY + ngj, 1.0));
    v[j]  = xtemp / (xtemp + ytemp);

    if(j != 0){
      w[j] = v[j] * ((1 - v[j-1]) * w[j - 1]) / v[j-1];
    }else{
      w[j] = v[j];
    }
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

void hyper_accelerate_SLI_mv_MRK(arma::vec y,
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

  double tval = 0.0;
  arma::vec beta_acc(beta0.n_elem, arma::fill::zeros);
  arma::mat Sb_acc(beta0.n_elem, beta0.n_elem, arma::fill::zeros);
  int k = 0;

  for(arma::uword j = 0; j < beta.n_rows; j++){
    if(arma::accu(clust == j) != 0){
      k += 1;
      beta_acc += arma::trans(beta.row(j));
      tval += 1 / sigma2(j);
    }
  }

  for(arma::uword j = 0; j < beta.n_rows; j++){
    if(arma::accu(clust == j) != 0){
      Sb_acc += (arma::trans(beta.row(j)) - beta_acc / k) * (beta.row(j) - arma::trans(beta_acc / k));
    }
  }

  // sampling hyperparameters
  double kn = k1 + k;
  arma::vec betan = ((beta1 * k1) + beta_acc)/kn;
  double sbn = sb1 + k;
  arma::mat Sbn = Sb1 + Sb_acc + ((k1 * k) / kn) *
    (beta_acc / k - beta1) * arma::trans(beta_acc / k - beta1);
  double tau1n = tau1 + k * a0;
  double tau2n = tau2 + tval;

  Sb0 = arma::inv(arma::wishrnd(arma::inv(Sbn), sbn));
  beta0 = arma::mvnrnd(betan, Sb0 / kn);
  b0 = arma::randg(arma::distr_param(tau1n, 1 / tau2n));

}

/*==================================================================================
 * grow parameters - MULTIVARIATE slice sampler - MRK
 * growing up the parameter vectors
 * till reaching the condition sum(w) > u_i, for all i
 *
 * args:
 * - mu:        matrix, each row a mean
 * - s2:        cube, each slice a covariance matrix
 * - v:         vector of stick break components
 * - w:         vector of stick weights
 * - u:         vector of uniform values
 * - m0:        vector of location's prior distribution
 * - k0:        vector of location's variance tuning parameter, one for each dimension
 * - a0:        vector of parameters for the scale component, one for each dimension
 * - b0:        vector of parameters for the scale component, one for each dimension
 * - mass:      mass parameter
 * - n:         number of observations
 * - sigma_PY:  discount parameter
 *
 * Void function
 ==================================================================================*/

void grow_param_SLI_PY_mv_MRK(arma::mat &beta,
                              arma::vec &sigma2,
                              arma::vec &v,
                              arma::vec &w,
                              arma::vec u,
                              arma::vec beta0,
                              arma::mat Sb0,
                              double a0,
                              double b0,
                              double mass,
                              int n,
                              double sigma_PY){

  double xtemp, ytemp;
  double w_sum = arma::accu(w);
  int k_old = beta.n_rows;
  int k = w.n_elem;
  int k_new;

  while((sum((1 - u) < w_sum) < n)){

    k = w.n_elem;
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

  k_new = w.n_elem;

  beta.resize(k_new, beta.n_cols);
  sigma2.resize(k_new);

  for(arma::uword j = k_old; j < k_new; j++){
    sigma2(j) = 1.0 / arma::randg(arma::distr_param(a0, 1.0 / (b0)));
    beta.row(j) = arma::trans(arma::mvnrnd(beta0, Sb0));
  }
}


/*==================================================================================
 * grow parameters - MULTIVARIATE slice sampler - MRK
 * growing up the parameter vectors
 * till reaching the condition sum(w) > u_i, for all i
 *
 * args:
 * - mu:        matrix, each row a mean
 * - s2:        cube, each slice a covariance matrix
 * - v:         vector of stick break components
 * - w:         vector of stick weights
 * - u:         vector of uniform values
 * - m0:        vector of location's prior distribution
 * - k0:        vector of location's variance tuning parameter, one for each dimension
 * - a0:        vector of parameters for the scale component, one for each dimension
 * - b0:        vector of parameters for the scale component, one for each dimension
 * - mass:      mass parameter
 * - n:         number of observations
 * - sigma_PY:  discount parameter
 *
 * Void function
 ==================================================================================*/

void grow_param_indep_SLI_PY_mv_MRK(arma::mat &beta,
                                    arma::vec &sigma2,
                                    arma::vec &v,
                                    arma::vec &w,
                                    arma::vec &xi,
                                    arma::vec u,
                                    arma::vec beta0,
                                    arma::mat Sb0,
                                    double a0,
                                    double b0,
                                    double mass,
                                    int n,
                                    double sigma_PY,
                                    double param_seq_one,
                                    double param_seq_two){

  double xtemp, ytemp;
  double xi_sum = arma::accu(xi);
  int k_old = beta.n_rows;
  int k = w.n_elem;
  int k_new;

  while((sum((1 - u) <= xi_sum) < n)){

    k = w.n_elem;
    v.resize(k+1);
    w.resize(k+1);
    xi.resize(k+1);

    xtemp = arma::randg(arma::distr_param(1.0 - sigma_PY, 1.0));
    ytemp = arma::randg(arma::distr_param(mass + (k + 1) * sigma_PY, 1.0));
    v[k]  = xtemp / (xtemp + ytemp);

    if(k == 0){
      w[k] = v[k];
    }else{
      w[k] = v[k] * ((1 - v[k-1]) * w[k - 1]) / v[k-1];
    }

    xi[k] = xi[k - 1] * (param_seq_one + k * param_seq_two) / (param_seq_one + 1 + k * param_seq_two);
    xi_sum += xi[k];
  }

  k_new = w.n_elem;

  beta.resize(k_new, beta.n_cols);
  sigma2.resize(k_new);

  for(arma::uword j = k_old; j < k_new; j++){
    sigma2(j) = 1.0 / arma::randg(arma::distr_param(a0, 1.0 / (b0)));
    beta.row(j) = arma::trans(arma::mvnrnd(beta0, Sb0));
  }
}

/*==================================================================================
 * Update clusters - MULTIVARIATE slice sampler - MRK
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

void update_cluster_SLI_mv_MRK(arma::vec y,
                               arma::mat covs,
                               arma::mat beta,
                               arma::vec sigma2,
                               arma::vec &clust,
                               arma::vec w,
                               arma::vec u){
  int n = covs.n_rows;
  int k = beta.n_rows;

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
    } else {
      probs.resize(index_use.n_elem);
      for(arma::uword j = 0; j < index_use.n_elem; j++){
        probs(j) = - log(2 * M_PI * sigma2(index_use(j))) / 2 -
          (pow(y(i) - arma::dot(covs.row(i), beta.row(index_use(j))), 2) /
            (2 * sigma2(index_use(j))));
      }

      sampled = rintnunif_log(probs);
      clust[i] = index_use[sampled];
    }
  }
}

/*==================================================================================
 * Update clusters - MULTIVARIATE slice sampler - MRK - INDEP
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

void update_cluster_indep_SLI_mv_MRK(arma::vec y,
                                     arma::mat covs,
                                     arma::mat beta,
                                     arma::vec sigma2,
                                     arma::vec &clust,
                                     arma::vec w,
                                     arma::vec xi,
                                     arma::vec u){
  int n = covs.n_rows;
  int k = beta.n_rows;

  arma::uvec index_use;
  arma::uvec index = arma::regspace<arma::uvec>(0, k - 1);
  arma::vec probs;
  int siz;
  int sampled;

  for(arma::uword i = 0; i < n; i++){
    siz = 0;
    index_use.resize(1);
    for(arma::uword r = 0; r < k; r++){
      if(xi[r] > u[i]){
        siz++;
        index_use.resize(siz);
        index_use[siz - 1] = index[r];
      }
    }

    if(index_use.n_elem == 1){
      clust[i] = index_use[0];
    } else {
      probs.resize(index_use.n_elem);
      for(arma::uword j = 0; j < index_use.n_elem; j++){
        probs(j) = log(w(index_use[j])) - log(xi(index_use[j])) - log(2 * M_PI * sigma2(index_use(j))) / 2 -
          (pow(y(i) - arma::dot(covs.row(i), beta.row(index_use(j))), 2) /
            (2 * sigma2(index_use(j))));
      }

      sampled = rintnunif_log(probs);
      clust[i] = index_use[sampled];
    }
  }
}

/*----------------------------------------------------------------------
 *
 * MIXTURE OF REGRESSION KERNELS
 * LOCATION KERNEL
 * SLI functions
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

void accelerate_SLI_mv_MRK_L(arma::vec y,
                             arma::mat covs,
                             arma::mat &beta,
                             double &sigma2,
                             arma::vec &v,
                             arma::vec &w,
                             arma::vec clust,
                             arma::vec beta0,
                             arma::mat Sb0,
                             double a0,
                             double b0,
                             double mass,
                             double sigma_PY){

  arma::mat cdata;
  arma::mat tdata;
  arma::vec ty;

  double an, bn, xtemp, ytemp;
  arma::mat tSb;
  arma::vec tbeta0;
  double accu_cdata = 0.0;

  int nj, ngj;

  // loop over the clusters
  for(arma::uword j = 0; j < beta.n_rows; j++){

    // INITIALIZE CLUSTER RELATED QUANTITIES
    nj  = arma::accu(clust == j);
    ngj = arma::accu(clust > j);
    tdata = covs.rows(arma::find(clust == j));
    ty = y.elem(arma::find(clust == j));

    // UPDATE THE LOCATIONS MARGINALLY
    if(tdata.n_rows > 0){
      cdata = ty - (tdata * beta.row(j).t());
      accu_cdata += arma::accu(pow(cdata, 2));

      // update the coefficients
      tSb = arma::inv(arma::inv(Sb0) + arma::trans(tdata) * tdata / sigma2);
      tbeta0 = tSb * (arma::inv(Sb0) * beta0 + (arma::trans(tdata) * ty) / sigma2);
      beta.row(j) = arma::trans(arma::mvnrnd(tbeta0, tSb));

    } else {

      beta.row(j) = arma::trans(arma::mvnrnd(beta0, Sb0));
    }

    // update the variance
    an = a0 + y.n_elem/2;
    bn = b0 + accu_cdata / 2;
    sigma2 = 1.0 / arma::randg(arma::distr_param(an, 1.0 / (bn)));


    // SAMPLE THE WEIGHTS OF STICK-BREAKING REPRESENTATION
    xtemp = arma::randg(arma::distr_param(1.0 - sigma_PY + nj, 1.0));
    ytemp = arma::randg(arma::distr_param(mass + (j + 1) * sigma_PY + ngj, 1.0));
    v[j]  = xtemp / (xtemp + ytemp);

    if(j != 0){
      w[j] = v[j] * ((1 - v[j-1]) * w[j - 1]) / v[j-1];
    }else{
      w[j] = v[j];
    }
  }
}

/*==================================================================================
 * Hyper-accelerate - MULTIVARIATE importance conditional sampler - MRK
 *
 * args:
 * - mu:      vector of location component
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

void hyper_accelerate_SLI_mv_MRK_L(arma::vec y,
                                   arma::mat covs,
                                   arma::vec clust,
                                   arma::mat beta,
                                   arma::vec &beta0,
                                   arma::mat &Sb0,
                                   arma::vec beta1,
                                   double k1,
                                   double sb1,
                                   arma::mat Sb1){

  arma::vec beta_acc(beta0.n_elem, arma::fill::zeros);
  arma::mat Sb_acc(beta0.n_elem, beta0.n_elem, arma::fill::zeros);
  int k = 0;

  for(arma::uword j = 0; j < beta.n_rows; j++){
    if(arma::accu(clust == j) != 0){
      k += 1;
      beta_acc += arma::trans(beta.row(j));
    }
  }

  for(arma::uword j = 0; j < beta.n_rows; j++){
    if(arma::accu(clust == j) != 0){
      Sb_acc += (arma::trans(beta.row(j)) - beta_acc / k) * (beta.row(j) - arma::trans(beta_acc / k));
    }
  }

  // sampling hyperparameters
  double kn = k1 + k;
  arma::vec betan = ((beta1 * k1) + beta_acc)/kn;
  double sbn = sb1 + k;
  arma::mat Sbn = Sb1 + Sb_acc + ((k1 * k) / kn) *
    (beta_acc / k - beta1) * arma::trans(beta_acc / k - beta1);

  Sb0 = arma::inv(arma::wishrnd(arma::inv(Sbn), sbn));
  beta0 = arma::mvnrnd(betan, Sb0 / kn);
}

/*==================================================================================
 * grow parameters - MULTIVARIATE slice sampler - MRK
 * growing up the parameter vectors
 * till reaching the condition sum(w) > u_i, for all i
 *
 * args:
 * - mu:        matrix, each row a mean
 * - v:         vector of stick break components
 * - w:         vector of stick weights
 * - u:         vector of uniform values
 * - m0:        vector of location's prior distribution
 * - k0:        vector of location's variance tuning parameter, one for each dimension
 * - mass:      mass parameter
 * - n:         number of observations
 * - sigma_PY:  discount parameter
 *
 * Void function
 ==================================================================================*/

void grow_param_SLI_PY_mv_MRK_L(arma::mat &beta,
                                arma::vec &v,
                                arma::vec &w,
                                arma::vec u,
                                arma::vec beta0,
                                arma::mat Sb0,
                                double mass,
                                int n,
                                double sigma_PY){

  double xtemp, ytemp;
  double w_sum = arma::accu(w);
  int k_old = beta.n_rows;
  int k = w.n_elem;
  int k_new;

  while((sum((1 - u) < w_sum) < n)){

    k = w.n_elem;
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

  k_new = w.n_elem;
  beta.resize(k_new, beta.n_cols);

  for(arma::uword j = k_old; j < k_new; j++){
    beta.row(j) = arma::trans(arma::mvnrnd(beta0, Sb0));
  }
}


/*==================================================================================
 * grow parameters - MULTIVARIATE slice sampler - MRK
 * growing up the parameter vectors
 * till reaching the condition sum(w) > u_i, for all i
 *
 * args:
 * - mu:        matrix, each row a mean
 * - s2:        cube, each slice a covariance matrix
 * - v:         vector of stick break components
 * - w:         vector of stick weights
 * - u:         vector of uniform values
 * - m0:        vector of location's prior distribution
 * - k0:        vector of location's variance tuning parameter, one for each dimension
 * - a0:        vector of parameters for the scale component, one for each dimension
 * - b0:        vector of parameters for the scale component, one for each dimension
 * - mass:      mass parameter
 * - n:         number of observations
 * - sigma_PY:  discount parameter
 *
 * Void function
 ==================================================================================*/

void grow_param_indep_SLI_PY_mv_MRK_L(arma::mat &beta,
                                      arma::vec &v,
                                      arma::vec &w,
                                      arma::vec &xi,
                                      arma::vec u,
                                      arma::vec beta0,
                                      arma::mat Sb0,
                                      double mass,
                                      int n,
                                      double sigma_PY,
                                      double param_seq_one,
                                      double param_seq_two){

  double xtemp, ytemp;
  double xi_sum = arma::accu(xi);
  int k_old = beta.n_rows;
  int k = w.n_elem;
  int k_new;

  while((sum((1 - u) <= xi_sum) < n)){

    k = w.n_elem;
    v.resize(k+1);
    w.resize(k+1);
    xi.resize(k+1);

    xtemp = arma::randg(arma::distr_param(1.0 - sigma_PY, 1.0));
    ytemp = arma::randg(arma::distr_param(mass + (k + 1) * sigma_PY, 1.0));
    v[k]  = xtemp / (xtemp + ytemp);

    if(k == 0){
      w[k] = v[k];
    }else{
      w[k] = v[k] * ((1 - v[k-1]) * w[k - 1]) / v[k-1];
    }

    xi[k] = xi[k - 1] * (param_seq_one + k * param_seq_two) / (param_seq_one + 1 + k * param_seq_two);
    xi_sum += xi[k];
  }

  k_new = w.n_elem;
  beta.resize(k_new, beta.n_cols);

  for(arma::uword j = k_old; j < k_new; j++){
    beta.row(j) = arma::trans(arma::mvnrnd(beta0, Sb0));
  }
}

/*==================================================================================
 * Update clusters - MULTIVARIATE slice sampler - MRK
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

void update_cluster_SLI_mv_MRK_L(arma::vec y,
                                 arma::mat covs,
                                 arma::mat beta,
                                 double sigma2,
                                 arma::vec &clust,
                                 arma::vec w,
                                 arma::vec u){
  int n = covs.n_rows;
  int k = beta.n_rows;

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
    } else {
      probs.resize(index_use.n_elem);
      for(arma::uword j = 0; j < index_use.n_elem; j++){
        probs(j) = - log(2 * M_PI * sigma2) / 2 -
          (pow(y(i) - arma::dot(covs.row(i), beta.row(index_use(j))), 2) /
            (2 * sigma2));
      }

      sampled = rintnunif_log(probs);
      clust[i] = index_use[sampled];
    }
  }
}

/*==================================================================================
 * Update clusters - MULTIVARIATE slice sampler - MRK - INDEP
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

void update_cluster_indep_SLI_mv_MRK_L(arma::vec y,
                                       arma::mat covs,
                                       arma::mat beta,
                                       double sigma2,
                                       arma::vec &clust,
                                       arma::vec w,
                                       arma::vec xi,
                                       arma::vec u){
  int n = covs.n_rows;
  int k = beta.n_rows;

  arma::uvec index_use;
  arma::uvec index = arma::regspace<arma::uvec>(0, k - 1);
  arma::vec probs;
  int siz;
  int sampled;

  for(arma::uword i = 0; i < n; i++){
    siz = 0;
    index_use.resize(1);
    for(arma::uword r = 0; r < k; r++){
      if(xi[r] > u[i]){
        siz++;
        index_use.resize(siz);
        index_use[siz - 1] = index[r];
      }
    }

    if(index_use.n_elem == 1){
      clust[i] = index_use[0];
    } else {
      probs.resize(index_use.n_elem);
      for(arma::uword j = 0; j < index_use.n_elem; j++){
        probs(j) = log(w(index_use[j])) - log(xi(index_use[j])) - log(2 * M_PI * sigma2) / 2 -
          (pow(y(i) - arma::dot(covs.row(i), beta.row(index_use(j))), 2) /
            (2 * sigma2));
      }

      sampled = rintnunif_log(probs);
      clust[i] = index_use[sampled];
    }
  }
}
