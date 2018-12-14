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
* Accelerate - UNIVARIATE conditional Polya Urn Scheme
* acceleration step for reshuffling the parameters
* given an allocation
*
* args:
* - data:  vector of observation
* - group: vector, each element correspond to a group
* - zeta:  vector, each element correpsond to an urn
* - mu:    vector of location component
* - s2:    vector of scale component
* - clust: vector of allocation
* - m0:    mean of location's prior distribution (scalar)
* - k0:    parameter of NIG on scale of location distribution (scalar)
* - a0:    shape of prior Gamma distribution on the scale component (scalar)
* - b0:    rate of prior Gamma distribution on the scale component (scalar)
* - ngr:   number of groups
*
* Void function
* ==================================================================================
*/

void accelerate_DDP(arma::vec data,
                     arma::vec group,
                     arma::vec zeta,
                     arma::field<arma::vec> &mu,
                     arma::field<arma::vec> &s2,
                     arma::vec clust,
                     double m0,
                     double k0,
                     double a0,
                     double b0,
                     int ngr){

  // for any urn
  for(arma::uword g = 0; g <= ngr; g++){

    // if the urn is not empty
    if(any(zeta == g)){

      // for any block in the urn g
      for(arma::uword j = 0; j < mu(g).n_elem; j++){

        // check if the block is not empty
        if(any(zeta == g && clust == j)){

          // update the corresponding parameters
          int nj = sum(zeta == g && clust == j);
          arma::vec tdata = data.elem(arma::find(zeta == g && clust == j));

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
}

/*==================================================================================
* Clean parameter - UNIVARIATE conditional Polya Urn Scheme
* discard the middle not used values for the clusters and
* update the correspondent parameters.
*
* args:
* - mu:      vector, each element a mean
* - s2:      vector, each element a variance
* - clust:   vector, each (integer) value is the cluster of the corresp. obs
* - group:   vector, each element correspond to a group
* - zeta:    vector, each element correpsond to an urn
* - ngr:     number of groups
*
* Void function.
* ==================================================================================
*/

void para_clean_DDP(arma::field<arma::vec> &mu,
                     arma::field<arma::vec> &s2,
                     arma::vec &clust,
                     arma::vec group,
                     arma::vec zeta,
                     int ngr) {

  // loop over the non-empty urn
  for(arma::uword g = 0; g <= ngr; g++){
    if(any(zeta == g)){

      int k = mu(g).n_elem;

      // for all the used parameters
      for(arma::uword i = 0; i < k; i++){

        // if a cluster is empty
        if(!any(zeta == g && clust == i)){

          // find the last full cluster, then swap
          for(arma::uword j = k; j > i; j--){
            if(any(zeta == g && clust == j)){

              // SWAPPING!!
              clust.elem(arma::find(zeta == g && clust == j)).fill(i);

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
      if(any(zeta == g)){
        int m_ind = max(clust.elem(arma::find(zeta == g))) + 1;
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
* Simulate finite distribution - UNIVARIATE conditional Polya Urn Scheme
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
* - ngr:         number of groups
*
* Void function.
* ==================================================================================
*/

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
    mutemp(g) = arma::randn(k) % sqrt(k0 * s2temp(g)) + m0;
  }
}

/*==================================================================================
* Update clusters - UNIVARIATE conditional Polya Urn Scheme
*
* args:
* - data:        vector of observation
* - group:       vector, each element correspond to a group
* - zeta:        vector, each element correpsond to an urn
* - mujoin:      mean values for each component
* - s2join:      mean values for each component
* - probjoin:    mean values for each component
* - clust:       vector of allocation
* - max_val:     maximum old value
* - iter:        current iteration
* - ngr:         number of groups
* - new_val      counting the new values
*
* void function
* ==================================================================================
*/

void clust_update_DDP(arma::vec data,
                      arma::vec group,
                      arma::vec &zeta,
                      arma::field<arma::vec> mujoin_complete,
                      arma::field<arma::vec> s2join_complete,
                      arma::field<arma::vec> probjoin_complete,
                      arma::vec &clust,
                      arma::vec max_val,
                      int iter,
                      int ngr){
  arma::vec probs_upd;
  int index;


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
    }

    // sample the allocation, if new update new_val
    clust(i) = rintnunif(probs_upd);
    if(clust(i) >= max_val(index)){
      zeta(i) = 0;
      clust(i) = clust(i) - max_val(index);
    } else {
      zeta(i) = group(i);
    }
  }
}

/*==================================================================================
* Update Urns weights - UNIVARIATE conditional Polya Urn Scheme
*
* args:
* - w:     vector of weights
* - mass:  parameter mass of the DP
* - wei:   weight parameter for the mixture of gamma process
*
* void function
* ==================================================================================
*/

void update_w_DDP(arma::vec &w,
                  double mass,
                  double wei,
                  int n_approx_unif,
                  arma::vec zeta,
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
  double cost_G = exp((ngr * std::lgamma(mass * wei)) - (std::lgamma(ngr * mass * wei)));
  for(arma::uword g = 0; g < ngr; g++){
    temp_exp_up(g)   = mass * wei - 1 + arma::accu(group == g + 1 && zeta == g + 1);
    temp_exp_down(g) = mass * wei + 1 - arma::accu(group == g + 1 && zeta != g + 1);
  }

  // loop over each sampled value from the uniform
  for(arma::uword j = 0; j < n_approx_unif; j++){

    // compute the terms
    temp = 0.0;
    for(arma::uword g = 0; g < ngr; g++){
      temp += std::pow(uvals.row(j)(g), temp_exp_up(g)) / std::pow(uvals.row(j)(g), temp_exp_down(g));
    }

    // update j-th element of the probs vector
    imp_probs(j) = cost_G * temp * std::pow(1 + arma::accu(uvals.row(j) / (1 - uvals.row(j))), - ngr * mass * wei);
  }

  int index = rintnunif(imp_probs);
  w = arma::trans(uvals.row(index));
}
