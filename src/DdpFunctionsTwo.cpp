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
 *==================================================================================
 */

#include "RcppArmadillo.h"
#include "Distributions.h"
// [[Rcpp::depends("RcppArmadillo")]]

/*==================================================================================
 * Acceleration step - GM-DDP
 *
 * args:
 * - data:       vector of observation
 * - group:      vector, each element correspond to a group
 * - group_log:  vector, each element correpsond to an urn
 * - col_log:    vector, each element correpsond to a common urn for the first variable
 * - row_log:    vector, each element correpsond to a common urn for the second variable
 * - mu:         field of location component
 * - s2:         field of scale component
 * - muc:        field of location component for the common urn of the first variable
 * - s2c:        field of scale component or the common urn of the first variable
 * - mur:        field of location component for the common urn of the second variable
 * - s2r:        field of scale component or the common urn of the second variable
 * - clust:      vector of allocation
 * - m0:         mean of location's prior distribution (scalar)
 * - k0:         parameter of NIG on scale of location distribution (scalar)
 * - a0:         shape of prior Gamma distribution on the scale component (scalar)
 * - b0:         rate of prior Gamma distribution on the scale component (scalar)
 * - ngr:        number of groups
 * - ngrc:       number of groups of the first variable
 * - ngrr:       number of groups of the second variable
 *
 * Void function
 * ==================================================================================
 */

void accelerate_DDP2(arma::vec data,
                     arma::vec group,
                     arma::vec group_log,
                     arma::vec col_log,
                     arma::vec row_log,
                     arma::field<arma::vec> &mu,
                     arma::field<arma::vec> &s2,
                     arma::field<arma::vec> &muc,
                     arma::field<arma::vec> &s2c,
                     arma::field<arma::vec> &mur,
                     arma::field<arma::vec> &s2r,
                     arma::vec clust,
                     double m0,
                     double k0,
                     double a0,
                     double b0,
                     int ngr,
                     int ngrc,
                     int ngrr){

  // SPECIFIC PROCESSES
  // for any urn
  for(arma::uword g = 0; g < ngr; g++){

    // if the urn is not empty
    if(any(group_log == g)){

      // for any block in the urn g
      for(arma::uword j = 0; j < mu(g).n_elem; j++){

        // check if the block is not empty
        if(any(group_log == g && clust == j)){

          // update the corresponding parameters
          int nj = sum(group_log == g && clust == j);
          arma::vec tdata = data.elem(arma::find(group_log == g && clust == j));

          double kn = 1.0 / ( (1.0/k0) + nj);
          double mn = kn * ((m0/k0) + sum(tdata));
          double an = a0 + (nj / 2.0);
          double bn = b0 + (((pow(m0, 2)/ k0) + arma::accu(arma::pow(tdata, 2)) - (pow(mn, 2)/ kn)) / 2.0);

          s2(g)(j) = 1.0 / arma::randg(arma::distr_param(an, 1.0 / bn));
          mu(g)(j) = arma::randn() * sqrt(kn * s2(g)(j)) + mn;
        } else {

          // otherwise update from the prior
          s2(g)(j) = 1.0 / arma::randg(arma::distr_param(a0, 1.0 / b0));
          mu(g)(j) = arma::randn() * sqrt(k0 * s2(g)(j)) + m0;
        }

      }
    } else {

      // if is empty, resize to 0 the objects
      mu(g).resize(0);
      s2(g).resize(0);
    }
  }

  // FIRST VARIABLE PROCESSES
  // for any urn
  for(arma::uword g = 0; g < ngrc; g++){

    // if the urn is not empty
    if(any(group_log == ngr && col_log == g)){

      // for any block in the urn g
      for(arma::uword j = 0; j < muc(g).n_elem; j++){

        // check if the block is not empty
        if(any(group_log == ngr && col_log == g && clust == j)){

          // update the corresponding parameters
          int nj = sum(group_log == ngr && col_log == g && clust == j);
          arma::vec tdata = data.elem(arma::find(group_log == ngr && col_log == g && clust == j));

          double kn = 1.0 / ( (1.0/k0) + nj);
          double mn = kn * ((m0/k0) + sum(tdata));
          double an = a0 + (nj / 2.0);
          double bn = b0 + (((pow(m0, 2)/ k0) + arma::accu(arma::pow(tdata, 2)) - (pow(mn, 2)/ kn)) / 2.0);

          s2c(g)(j) = 1.0 / arma::randg(arma::distr_param(an, 1.0 / bn));
          muc(g)(j) = arma::randn() * sqrt(kn * s2c(g)(j)) + mn;
        } else {

          // otherwise update from the prior
          s2c(g)(j) = 1.0 / arma::randg(arma::distr_param(a0, 1.0 / b0));
          muc(g)(j) = arma::randn() * sqrt(k0 * s2c(g)(j)) + m0;
        }

      }
    } else {

      // if is empty, resize to 0 the objects
      muc(g).resize(0);
      s2c(g).resize(0);
    }
  }

  // SECOND VARIABLE PROCESSES
  // for any urn
  for(arma::uword g = 0; g < ngrr; g++){

    // if the urn is not empty
    if(any(group_log == ngr && row_log == g)){

      // for any block in the urn g
      for(arma::uword j = 0; j < mur(g).n_elem; j++){

        // check if the block is not empty
        if(any(group_log == ngr && row_log == g && clust == j)){

          // update the corresponding parameters
          int nj = sum(group_log == ngr && row_log == g && clust == j);
          arma::vec tdata = data.elem(arma::find(group_log == ngr && row_log == g && clust == j));

          double kn = 1.0 / ( (1.0/k0) + nj);
          double mn = kn * ((m0/k0) + sum(tdata));
          double an = a0 + (nj / 2.0);
          double bn = b0 + (((pow(m0, 2)/ k0) + arma::accu(arma::pow(tdata, 2)) - (pow(mn, 2)/ kn)) / 2.0);

          s2r(g)(j) = 1.0 / arma::randg(arma::distr_param(an, 1.0 / bn));
          mur(g)(j) = arma::randn() * sqrt(kn * s2r(g)(j)) + mn;
        } else {

          // otherwise update from the prior
          s2r(g)(j) = 1.0 / arma::randg(arma::distr_param(a0, 1.0 / b0));
          mur(g)(j) = arma::randn() * sqrt(k0 * s2r(g)(j)) + m0;
        }

      }
    } else {

      // if is empty, resize to 0 the objects
      mur(g).resize(0);
      s2r(g).resize(0);
    }
  }
}

/*==================================================================================
 * Clean parameter - GM-DDP
 *
 * args:
 * - mu:      field, locations of the specific processes
 * - s2:      field, dispersion of the specific processes
 * - muc:     field, locations of the first variable processes
 * - s2c:     field, dispersion of the first variable processes
 * - mur:     field, locations of the second variable processes
 * - s2r:     field, dispersion of the second variable processes
 * - clust:   vector, each (integer) value is the cluster of the corresp. obs
 * - group:   vector, each element correspond to a group
 * - group_log:    vector, each element correpsond to an urn
 * - col_log:    vector, each element correpsond to a common urn for the first variable
 * - row_log:    vector, each element correpsond to a common urn for the second variable
 * - ngr:     number of groups
 * - ngrc:       number of groups of the first variable
 * - ngrr:       number of groups of the second variable
 *
 * Void function.
 * ==================================================================================
 */

void para_clean_DDP2(arma::field<arma::vec> &mu,
                     arma::field<arma::vec> &s2,
                     arma::field<arma::vec> &muc,
                     arma::field<arma::vec> &s2c,
                     arma::field<arma::vec> &mur,
                     arma::field<arma::vec> &s2r,
                     arma::vec &clust,
                     arma::vec group,
                     arma::vec group_log,
                     arma::vec col_log,
                     arma::vec row_log,
                     int ngr,
                     int ngrc,
                     int ngrr) {

  // SPECIFIC PROCESSES
  // loop over the non-empty urn
  for(arma::uword g = 0; g < ngr; g++){
    if(any(group_log == g)){

      int k = mu(g).n_elem;

      // for all the used parameters
      for(arma::uword i = 0; i < k; i++){

        // if a cluster is empty
        if(!any(group_log == g && clust == i)){

          // find the last full cluster, then swap
          for(arma::uword j = k; j > i; j--){
            if(any(group_log == g && clust == j)){

              // swap the clust values,
              // location and scale parameters
              clust.elem(arma::find(group_log == g && clust == j)).fill(i);

              double tmu = mu(g)(i);
              mu(g)(i) = mu(g)(j);
              mu(g)(j) = tmu;

              double ts2 = s2(g)(i);
              s2(g)(i) = s2(g)(j);
              s2(g)(j) = ts2;

              break;
            }
          }
        }
      }

      // reduce dimensions to the maximum needed
      if(any(group_log == g)){
        int m_ind = max(clust.elem(arma::find(group_log == g))) + 1;
        mu(g).resize(m_ind);
        s2(g).resize(m_ind);
      }
    } else {

      // if the urn is empty, resize to 0
      mu(g).resize(0);
      s2(g).resize(0);
    }
  }

  // FIRST VARIABLE PROCESSES
  // loop over the non-empty urn
  for(arma::uword g = 0; g < ngrc; g++){
    if(any(group_log == ngr && col_log == g)){

      int k = muc(g).n_elem;

      // for all the used parameters
      for(arma::uword i = 0; i < k; i++){

        // if a cluster is empty
        if(!any(group_log == ngr && col_log == g && clust == i)){

          // find the last full cluster, then swap
          for(arma::uword j = k; j > i; j--){
            if(any(group_log == ngr && col_log == g && clust == j)){

              // swap the clust values,
              // location and scale parameters
              clust.elem(arma::find(group_log == ngr && col_log == g && clust == j)).fill(i);

              double tmuc = muc(g)(i);
              muc(g)(i) = muc(g)(j);
              muc(g)(j) = tmuc;

              double ts2c = s2c(g)(i);
              s2c(g)(i) = s2c(g)(j);
              s2c(g)(j) = ts2c;

              break;
            }
          }
        }
      }

      // reduce dimensions
      if(any(group_log == ngr && col_log == g)){
        int m_ind = max(clust.elem(arma::find(group_log == ngr && col_log == g))) + 1;
        muc(g).resize(m_ind);
        s2c(g).resize(m_ind);
      }
    } else {

      // if the urn is empty, resize to 0
      muc(g).resize(0);
      s2c(g).resize(0);
    }
  }

  // SECOND VARIABLE PROCESSES
  // loop over the non-empty urn - ROW
  for(arma::uword g = 0; g < ngrr; g++){
    if(any(group_log == ngr && row_log == g)){

      int k = mur(g).n_elem;

      // for all the used parameters
      for(arma::uword i = 0; i < k; i++){

        // if a cluster is empty
        if(!any(group_log == ngr && row_log == g && clust == i)){

          // find the last full cluster, then swap
          for(arma::uword j = k; j > i; j--){
            if(any(group_log == ngr && row_log == g && clust == j)){

              // swap the clust values,
              // location and scale parameters
              clust.elem(arma::find(group_log == ngr && row_log == g && clust == j)).fill(i);

              double tmur = mur(g)(i);
              mur(g)(i) = mur(g)(j);
              mur(g)(j) = tmur;

              double ts2r = s2r(g)(i);
              s2r(g)(i) = s2r(g)(j);
              s2r(g)(j) = ts2r;

              break;
            }
          }
        }
      }

      // reduce dimensions
      if(any(group_log == ngr && row_log == g)){
        int m_ind = max(clust.elem(arma::find(group_log == ngr && row_log == g))) + 1;
        mur(g).resize(m_ind);
        s2r(g).resize(m_ind);
      }
    } else {

      // if the urn is empty, resize to 0
      mur(g).resize(0);
      s2r(g).resize(0);
    }
  }
}

/*==================================================================================
 * Simulate finite distribution - GM-DDP
 *
 * args:
 * - mutemp:      field, mean values for each temp component, specific processes
 * - s2temp:      field, variance values for temp each component, specific processes
 * - freqtemp:    field, frequency values for temp each component, specific processes
 * - mutemp:      field, mean values for each temp component, first variable processes
 * - s2temp:      field, variance values for temp each component, first variable processes
 * - freqtemp:    field, frequency values for temp each component, first variable processes
 * - mutemp:      field, mean values for each temp component, second variable processes
 * - s2temp:      field, variance values for temp each component, second variable processes
 * - freqtemp:    field, frequency values for temp each component, second variable processes
 * - mass:        mass of Dirichlet process
 * - wei_group:   prior weight for the specific processes
 * - wei_col:     prior weight for the first variable processes
 * - m0:          mean of location's prior distribution (scalar)
 * - k0:          parameter of NIG on scale of location distribution (scalar)
 * - a0:          shape of prior Gamma distribution on the scale component (scalar)
 * - b0:          rate of prior Gamma distribution on the scale component (scalar)
 * - napprox:     number of values to be sampled
 * - ngr:         number of groups
 * - ngrc:        number of groups of the first variable
 * - ngrr:        number of groups of the second variable
 *
 * Void function.
 * ==================================================================================
 */

void simu_trunc_DDP2(arma::field<arma::vec> &mutemp,
                     arma::field<arma::vec> &s2temp,
                     arma::field<arma::vec> &freqtemp,
                     arma::field<arma::vec> &mutempc,
                     arma::field<arma::vec> &s2tempc,
                     arma::field<arma::vec> &freqtempc,
                     arma::field<arma::vec> &mutempr,
                     arma::field<arma::vec> &s2tempr,
                     arma::field<arma::vec> &freqtempr,
                     double mass,
                     double wei_group,
                     double wei_col,
                     double m0,
                     double k0,
                     double a0,
                     double b0,
                     int napprox,
                     int ngr,
                     int ngrc,
                     int ngrr){
  int temp_cl;

  // SPECIFIC PROCESSES
  // loop over the different urns
  for(arma::uword g = 0; g < ngr; g++){

    // By independence of atoms and jumps, sample
    // the jumps, then the entire vector of atoms

    // resize the objects to size 1
    freqtemp(g).resize(1);

    // initialize the first element
    freqtemp(g).fill(1);
    int k = 1;

    // generate napprox values with ties
    for(arma::uword j = 1; j < napprox; j++){

      temp_cl = rintnunifw(freqtemp(g), mass * wei_group);
      if(temp_cl < (k - 1)){

        // if is an existing one, increase the freq
        freqtemp(g)(temp_cl) += 1;

      } else {

        // if is a new one, generate the new parameters
        freqtemp(g).resize(k + 1);
        freqtemp(g)(k) = 1;
        k += 1;

      }
    }

    // resize the objects and generate the atoms
    mutemp(g).resize(k);
    s2temp(g).resize(k);
    s2temp(g) = 1.0 / arma::randg(k, arma::distr_param(a0, 1.0 / b0));
    mutemp(g) = arma::randn(k) % sqrt(k0 * s2temp(g)) + m0;
  }

  // FIRST VARIABLE PROCESSES
  // loop over the different urns
  for(arma::uword g = 0; g < ngrc; g++){

    // By independence of atoms and jumps, sample
    // the jumps, then the entire vector of atoms

    // resize the objects to size 1
    freqtempc(g).resize(1);

    // initialize the first element
    freqtempc(g).fill(1);
    int k = 1;

    // generate napprox values with ties
    for(arma::uword j = 1; j < napprox; j++){

      temp_cl = rintnunifw(freqtempc(g), mass * (1 - wei_group) * wei_col);
      if(temp_cl < (k - 1)){

        // if is an existing one, increase the freq
        freqtempc(g)(temp_cl) += 1;

      } else {

        // if is a new one, generate the new parameters
        freqtempc(g).resize(k + 1);
        freqtempc(g)(k) = 1;
        k += 1;

      }
    }

    // resize the objects and generate the atoms
    mutempc(g).resize(k);
    s2tempc(g).resize(k);
    s2tempc(g) = 1.0 / arma::randg(k, arma::distr_param(a0, 1.0 / b0));
    mutempc(g) = arma::randn(k) % sqrt(k0 * s2tempc(g)) + m0;
  }

  // SECOND VARIABLE PROCESSES
  // loop over the different urns
  for(arma::uword g = 0; g < ngrr; g++){

    // By independence of atoms and jumps, sample
    // the jumps, then the entire vector of atoms

    // resize the objects to size 1
    freqtempr(g).resize(1);

    // initialize the first element
    freqtempr(g).fill(1);
    int k = 1;

    // generate napprox values with ties
    for(arma::uword j = 1; j < napprox; j++){

      temp_cl = rintnunifw(freqtempr(g), mass * (1 - wei_group) * (1 - wei_col));
      if(temp_cl < (k - 1)){

        // if is an existing one, increase the freq
        freqtempr(g)(temp_cl) += 1;

      } else {

        // if is a new one, generate the new parameters
        freqtempr(g).resize(k + 1);
        freqtempr(g)(k) = 1;
        k += 1;

      }
    }

    // resize the objects and generate the atoms
    mutempr(g).resize(k);
    s2tempr(g).resize(k);
    s2tempr(g) = 1.0 / arma::randg(k, arma::distr_param(a0, 1.0 / b0));
    mutempr(g) = arma::randn(k) % sqrt(k0 * s2tempr(g)) + m0;
  }
}

/*==================================================================================
 * Update clusters - GM-DDP
 *
 * args:
 * - data:          vector of observation
 * - group:         vector, each element correspond to a group
 * - col_group:     vector, each element correspond to a group for the first variable
 * - row_group:     vector, each element correspond to a group for the second variable
 * - group_log:     vector, each element correpsond to an urn
 * - col_log:       vector, each element correpsond to a common urn for the first variable
 * - row_log:       vector, each element correpsond to a common urn for the second variable
 * - mujoin:        field, mean values for each component
 * - s2join:        field, mean values for each component
 * - probjoin:      field, mean values for each component
 * - clust:         vector of allocation
 * - max_val_group: number of atoms from the specific processes
 * - max_val_col:   number of atoms from the first variable processes
 * - iter:          current iteration
 * - ngr:           number of groups
 * - ngrc:          number of groups of the first variable
 * - ngrr:          number of groups of the second variable
 *
 * void function
 * ==================================================================================
 */

void clust_update_DDP2(arma::vec data,
                       arma::vec group,
                       arma::vec col_group,
                       arma::vec row_group,
                       arma::vec &group_log,
                       arma::vec &col_log,
                       arma::vec &row_log,
                       arma::field<arma::vec> mujoin_complete,
                       arma::field<arma::vec> s2join_complete,
                       arma::field<arma::vec> probjoin_complete,
                       arma::vec &clust,
                       arma::vec max_val_group,
                       arma::vec max_val_col,
                       int ngr,
                       int ngrc,
                       int ngrr){
  arma::vec probs_upd;
  int index;
  int indexc;

  // for each observation
  for(arma::uword i = 0; i < data.n_elem; i++){
    index  = group(i);
    indexc = col_group(i);

    // resize the probs_upd object as the number of potential blocks
    int k = probjoin_complete(index).n_elem;
    probs_upd.resize(k);

    // compute the probability of each block
    for(arma::uword j = 0; j < k; j++){

      probs_upd(j) = probjoin_complete(index)(j) *
        arma::normpdf(data(i), mujoin_complete(index)(j),
                      sqrt(s2join_complete(index)(j)));

    }

    // sample the allocation, if new update new_val
    clust(i) = rintnunif(probs_upd);
    if(clust(i) < max_val_group(index)){

      group_log(i) = group(i);
      col_log(i) = ngrc;
      row_log(i) = ngrr;

    } else if((clust(i) >= max_val_group(index)) &&
      (clust(i) < (max_val_group(index) + max_val_col(indexc)))) {

      group_log(i) = ngr;
      col_log(i) = col_group(i);
      row_log(i) = ngrr;
      clust(i) = clust(i) - max_val_group(index);

    } else {

      group_log(i) = ngr;
      col_log(i) = ngrc;
      row_log(i) = row_group(i);
      clust(i) = clust(i) - max_val_group(index) - max_val_col(indexc);

    }
  }
}

/*==================================================================================
 * Update Urns weights - GM-DDP
 *
 * args:
 * - w:              vector of weights
 * - mass:           parameter mass of the DP
 * - wei_group:      prior weight for the specific processes
 * - n_approx_unif:  number of values to be sampled form the uniform
 * - group_log:      vector, each element correpsond to an urn
 * - clust:          vector of allocation
 * - group:          vector, each element correspond to a group
 * - ngr:            number of groups
 *
 * void function
 * ==================================================================================
 */

void update_w_DDP2(arma::vec &w,
                   double mass,
                   double wei_group,
                   int n_approx_unif,
                   arma::vec group_log,
                   arma::vec clust,
                   arma::vec group,
                   int ngr){

  // generate n_approx_unif values from an uniform distribution defined
  // on the hypercube with same dimension as w
  // and the vector for the probabilities
  arma::mat uvals(n_approx_unif, ngr, arma::fill::randu);
  arma::vec imp_probs(n_approx_unif);
  arma::vec temp_exp_up(ngr);
  arma::vec temp_exp_down(ngr);
  arma::vec temp_vals(ngr);
  double temp;

  // initialize the Beta costant
  // initialize the exponents (same for each simulated unif value)
  double cost_G = exp((ngr * std::lgamma(mass * wei_group)) - (std::lgamma(ngr * mass * wei_group)));
  for(arma::uword g = 0; g < ngr; g++){
    temp_exp_up(g)   = mass * wei_group - 1 + arma::accu(group == g && group_log == g);
    temp_exp_down(g) = mass * wei_group + 1 - arma::accu(group == g && group_log != g);
  }

  // loop over each sampled value from the uniform
  for(arma::uword j = 0; j < n_approx_unif; j++){

    // compute the terms
    temp = 0.0;
    for(arma::uword g = 0; g < ngr; g++){
      temp += std::pow(uvals.row(j)(g), temp_exp_up(g)) / std::pow(uvals.row(j)(g), temp_exp_down(g));
    }

    // update j-th element of the probs vector
    imp_probs(j) = cost_G * temp * std::pow(1 + arma::accu(uvals.row(j) / (1 - uvals.row(j))), - ngr * mass * wei_group);
  }

  int index = rintnunif(imp_probs);
  w = arma::trans(uvals.row(index));
}

/*==================================================================================
 * Update first variable weights - GM-DDP
 *
 * args:
 * - v:              vector of weights
 * - mass:           parameter mass of the DP
 * - wei_group:      prior weight for the specific processes
 * - wei_col:        prior weight for the first variable processes
 * - n_approx_unif:  number of values to be sampled form the uniform
 * - col_log:        vector, each element correpsond to a common urn for the first variable
 * - clust:          vector of allocation
 * - group:          vector, each element correspond to a group
 * - col_group:     vector, each element correspond to a group for the first variable
 * - ngr:            number of groups
 * - ngrc:           number of groups of the first variable
 *
 * void function
 * ==================================================================================
 */


void update_v_DDP2(arma::vec &v,
                   double mass,
                   double wei_group,
                   double wei_col,
                   int n_approx_unif,
                   arma::vec col_log,
                   arma::vec clust,
                   arma::vec group,
                   arma::vec col_group,
                   int ngr,
                   int ngrc){

  // generate n_approx_unif values from an uniform distribution defined
  // on the hypercube with same dimension as w
  // and the vector for the probabilities
  arma::mat uvals(n_approx_unif, ngrc, arma::fill::randu);
  arma::vec imp_probs(n_approx_unif);
  arma::vec temp_exp_up(ngrc);
  arma::vec temp_exp_down(ngrc);
  arma::vec temp_vals(ngrc);
  double temp;

  // initialize the Beta costant
  // initialize the exponents (same for each simulated unif value)
  double cost_G = exp((ngrc * std::lgamma(mass * (1 - wei_group) * wei_col)) -
                      (std::lgamma(ngrc * mass * (1 - wei_group) * wei_col)));
  for(arma::uword g = 0; g < ngrc; g++){
    temp_exp_up(g)   = mass * (1 - wei_group) * wei_col - 1 +
                        arma::accu(group == ngr && col_group == g && col_log == g);
    temp_exp_down(g) = mass * (1 - wei_group) * wei_col + 1 -
                        arma::accu(group == ngr && col_group == g && col_log != g);
  }

  // loop over each sampled value from the uniform
  for(arma::uword j = 0; j < n_approx_unif; j++){

    // compute the terms
    temp = 0.0;
    for(arma::uword g = 0; g < ngrc; g++){
      temp += std::pow(uvals.row(j)(g), temp_exp_up(g)) / std::pow(uvals.row(j)(g), temp_exp_down(g));
    }

    // update j-th element of the probs vector
    imp_probs(j) = cost_G * temp * std::pow(1 + arma::accu(uvals.row(j) /
                     (1 - uvals.row(j))), - ngrc * mass * (1 - wei_group) * wei_col);
  }

  int index = rintnunif(imp_probs);
  v = arma::trans(uvals.row(index));
}

