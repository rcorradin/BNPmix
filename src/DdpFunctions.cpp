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
 *==================================================================================*/

#include "RcppArmadillo.h"
#include "Distributions.h"
// [[Rcpp::depends("RcppArmadillo")]]

/*==================================================================================
 * Accelerate - UNIVARIATE importance conditional sampler for DDP model
 * acceleration step for reshuffling the parameters
 * given an allocation
 *
 * args:
 * - data:        vector of observation
 * - group:       vector, each element correspond to a group
 * - group_log:   vector, each element correpsond to group or the common one
 * - mu:          field of vector of location component
 * - s2:          field of vector of scale component
 * - clust:       vector of allocation
 * - m0:          mean of location's prior distribution (scalar)
 * - k0:          parameter of NIG on scale of location distribution (scalar)
 * - a0:          shape of prior Gamma distribution on the scale component (scalar)
 * - b0:          rate of prior Gamma distribution on the scale component (scalar)
 * - ngr:         number of groups
 *
 * Void function
 * ==================================================================================*/

void accelerate_DDP(arma::vec data,
                     arma::vec group,
                     arma::vec group_log,
                     arma::field<arma::vec> &mu,
                     arma::field<arma::vec> &s2,
                     arma::vec clust,
                     double m0,
                     double k0,
                     double a0,
                     double b0,
                     int ngr){
  double kn, mn, an, bn, data_m;
  arma::vec tdata;
  int nj;

  // for any urn
  for(arma::uword g = 0; g <= ngr; g++){

    // if the urn is not empty
    if(any(group_log == g)){

      // for any block in the urn g
      for(arma::uword j = 0; j < mu(g).n_elem; j++){

        // check if the block is not empty
        if(any(group_log == g && clust == j)){

          // update the corresponding parameters
          nj = sum(group_log == g && clust == j);
          tdata = data.elem(arma::find(group_log == g && clust == j));
          data_m = arma::accu(tdata) / nj;

          kn = (k0 + nj);
          mn = ((m0 * k0) + nj * data_m) / kn;
          an = a0 + (nj / 2.0);
          bn = b0 + (arma::accu(pow(tdata - data_m, 2)) + (nj * k0 * pow(data_m - m0, 2)) / (kn)) / 2;

          s2(g)(j) = 1.0 / arma::randg(arma::distr_param(an, 1.0 / bn));
          mu(g)(j) = arma::randn() * sqrt(s2(g)(j)/kn) + mn;
        } else {

          // otherwise update from the prior
          s2(g)(j) = 1.0 / arma::randg(arma::distr_param(a0, 1.0 / b0));
          mu(g)(j) = arma::randn() * sqrt(s2(g)(j) / k0) + m0;
        }

      }
    } else {

      // if is empty, resize to 0 the objects
      mu(g).resize(0);
      s2(g).resize(0);
    }
  }
}

/*==================================================================================
 * Clean parameter - UNIVARIATE importance conditional sampler for DDP models
 * discard the middle not used values for the clusters and
 * update the correspondent parameters.
 *
 * args:
 * - mu:          field of vector of location component
 * - s2:          field of vector of scale component
 * - clust:       vector, each (integer) value is the cluster of the corresp. obs
 * - group:       vector, each element correspond to a group
 * - group_log:   vector, each element correpsond to a group or the common one
 * - ngr:         number of groups
 *
 * Void function.
 * ==================================================================================*/

void para_clean_DDP(arma::field<arma::vec> &mu,
                     arma::field<arma::vec> &s2,
                     arma::vec &clust,
                     arma::vec group,
                     arma::vec group_log,
                     int ngr) {

  // loop over the non-empty urn
  for(arma::uword g = 0; g <= ngr; g++){
    if(any(group_log == g)){

      int k = mu(g).n_elem;

      // for all the used parameters
      for(arma::uword i = 0; i < k; i++){

        // if a cluster is empty
        if(!any(group_log == g && clust == i)){

          // find the last full cluster, then swap
          for(arma::uword j = k; j > i; j--){
            if(any(group_log == g && clust == j)){

              // SWAPPING!!
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

      // reduce dimensions
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
}

/*==================================================================================
 * Simulate finite distribution - importance conditional sampler for DDP models
 *
 * args:
 * - mutemp:      field of vector of mean values for each temp component
 * - s2temp:      field of vector of variance values for temp each component
 * - freqtemp:    field of vector of frequency values for temp each component
 * - mass:        mass of Dirichlet process
 * - wei:         vector of weights of the processes
 * - m0:          mean of location's prior distribution (scalar)
 * - k0:          parameter of NIG on scale of location distribution (scalar)
 * - a0:          shape of prior Gamma distribution on the scale component (scalar)
 * - b0:          rate of prior Gamma distribution on the scale component (scalar)
 * - napprox:     number of approximating values
 * - ngr:         number of groups
 *
 * Void function.
 * ==================================================================================*/

void simu_trunc_DDP(arma::field<arma::vec> &mutemp,
                    arma::field<arma::vec> &s2temp,
                    arma::field<arma::vec> &freqtemp,
                    double mass,
                    double wei,
                    double m0,
                    double k0,
                    double a0,
                    double b0,
                    int napprox,
                    int ngr){
  int temp_cl;

  // loop over the different urns
  for(arma::uword g = 0; g <= ngr; g++){

    // By independence of atoms and jumps, sample
    // the jumps, then the entire vector of atoms

    // resize the objects to size 1
    freqtemp(g).resize(1);

    // initialize the first element
    freqtemp(g).fill(1);
    int k = 1;

    // generate napprox values with ties
    for(arma::uword j = 1; j < napprox; j++){

      if(g != 0){
        temp_cl = rintnunifw(freqtemp(g), mass * wei);
      } else {
        temp_cl = rintnunifw(freqtemp(g), mass * (1 - wei));
      }

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
    mutemp(g) = arma::randn(k) % sqrt(s2temp(g) / k0) + m0;
  }
}

/*==================================================================================
 * Update clusters - UNIVARIATE importance conditional sampler for DDP models
 *
 * args:
 * - data:          vector of observation
 * - group:         vector, each element correspond to a group
 * - group_log:     vector, each element correpsond to a group or the common one
 * - mujoin:        field of vectors of mean values for each component
 * - s2join:        field of vectors of mean values for each component
 * - probjoin:      field of vectors of mean values for each component
 * - clust:         vector of allocation
 * - temp_proc_cum: cumulate values for processes
 * - max_val:       maximum old value
 * - iter:          current iteration
 * - ngr:           number of groups
 *
 * void function
 * ==================================================================================*/

void clust_update_DDP(arma::vec data,
                      arma::vec group,
                      arma::vec &group_log,
                      arma::field<arma::vec> mujoin_complete,
                      arma::field<arma::vec> s2join_complete,
                      arma::field<arma::vec> probjoin_complete,
                      arma::vec &clust,
                      arma::mat &temp_proc_cum,
                      arma::vec max_val,
                      int iter,
                      int ngr){
  arma::vec probs_upd;
  int index;
  temp_proc_cum.fill(0.0);

  // for each observation
  for(arma::uword i = 0; i < data.n_elem; i++){
    index = group(i) - 1;

    // resize the probs_upd object as the number of potential blocks
    int k = probjoin_complete(index).n_elem;
    probs_upd.resize(k);

    // compute the probability of each block
    for(arma::uword j = 0; j < k; j++){
      probs_upd(j) = probjoin_complete(index)(j) *
        arma::normpdf(data(i), mujoin_complete(index)(j),
                      sqrt(s2join_complete(index)(j)));

      if(j < max_val(index)){
        temp_proc_cum(i, 0) += probs_upd(j);
      } else {
        temp_proc_cum(i, 1) += probs_upd(j);
      }

    }

    // sample the allocation, if new update new_val
    clust(i) = rintnunif(probs_upd);
    if(clust(i) >= max_val(index)){
      group_log(i) = 0;
      clust(i) = clust(i) - max_val(index);
    } else {
      group_log(i) = group(i);
    }
  }
}

/*==================================================================================
* Update Urns weights - importance conditional sampler for DDP models
*
* args:
* - w:              vector of weights
* - mass:           parameter mass of the DP
* - wei:            weight of the processes
* - n_approx_unif:  number of values for the importance sampling step
* - group_log:      vector, each element a group or the common one
* - clust:          vector, clusters
* - group:          vector, each element a group
* - temp_proc_cum:  cumulated values for the processes for each observation
* - ngr:            number of groups
*
* void function
* ==================================================================================
*/

void update_w_DDP(arma::vec &w,
                  double mass,
                  double wei,
                  int n_approx_unif,
                  arma::vec group_log,
                  arma::vec clust,
                  arma::vec group,
                  arma::mat temp_proc_cum,
                  int ngr){

  // generate n_approx_unif values from an uniform distribution defined
  // on the hypercube with same dimension as w
  // and the vector for the probabilities
  int n = group.n_elem;
  arma::mat uvals(n_approx_unif, ngr, arma::fill::randu);
  // arma::mat uvals(49*49, 2);
  // arma::vec tv = arma::regspace( 0.02, 0.02, 1 );
  // int k = 0;
  // for(arma::uword i = 0; i < 49; i++){
  //   for(arma::uword j = 0; j < 49; j++){
  //     uvals(k,0) = tv(i);
  //     uvals(k,1) = tv(j);
  //     k++;
  //   }
  // }

  arma::vec imp_probs(uvals.n_rows);

  double tempval;

  // loop over each sampled value from the uniform
  // for(arma::uword j = 0; j < n_approx_unif; j++){
  for(arma::uword j = 0; j < uvals.n_rows; j++){

    tempval = 0.0;
    for(arma::uword g = 0; g < ngr; g++){

      tempval += log(std::pow(uvals.row(j)(g), mass * wei - 1) /
        std::pow(1 - uvals.row(j)(g), mass * wei + 1));

      for(arma::uword i = 0; i < n; i++){
        if(group(i) == g){
          tempval += log(uvals.row(j)(g) * temp_proc_cum(i,0) +
            (1 - uvals.row(j)(g)) * temp_proc_cum(i,1));
        }
      }

      tempval += (- ngr * mass * wei - mass * (1 - wei)) * log(1 + arma::accu(uvals.row(j) / (1 - uvals.row(j))));
    }


    // update j-th element of the probs vector
    imp_probs(j) = tempval;
  }

  int index = rintnunif_log(imp_probs);
  w = arma::trans(uvals.row(index));

  // // generate n_approx_unif values from an uniform distribution defined
  // // on the hypercube with same dimension as w
  // // and the vector for the probabilities
  // int n = group.n_elem;
  // double tempval;
  // arma::vec temp_w(w);
  // double const_w = 0.0;
  // arma::vec uvals(n_approx_unif);
  // arma::vec imp_probs(n_approx_unif);
  // // arma::vec uvals = arma::regspace( 0.02, 0.02, 0.98 );
  // // arma::vec imp_probs(uvals.n_elem);
  //
  // // loop over each sampled value from the uniform
  // for(arma::uword g = 0; g < ngr; g++){
  //
  //   uvals.randu();
  //   imp_probs.fill(0.0);
  //   temp_w = w;
  //   temp_w(g) = 0;
  //   const_w = 1 + arma::accu(temp_w / (1 - temp_w));
  //
  //   // sample the g element of w
  //   for(arma::uword j = 0; j < uvals.n_elem; j++){
  //
  //     tempval = 0.0;
  //     tempval += log(std::pow(uvals(j), mass * wei - 1) / std::pow(1 - uvals(j), mass * wei + 1));
  //
  //     for(arma::uword i = 0; i < n; i++){
  //       if(group(i) == g){
  //         tempval += log(uvals(j) * temp_proc_cum(i,0) +
  //           (1 - uvals(j)) * temp_proc_cum(i,1));
  //       }
  //     }
  //
  //     tempval += (- ngr * mass * wei - mass * (1 - wei)) * log(const_w + (uvals(j) / ( 1 - uvals(j))));
  //     imp_probs(j) = tempval;
  //   }
  //
  //   int index = rintnunif_log(imp_probs);
  //   w(g) = uvals(index);
  //
  // }
}

