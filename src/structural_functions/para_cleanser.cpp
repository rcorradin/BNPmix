/*
  Copyright (C) 2017 Riccardo Corradin

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

#include <distributions/gaussian.hpp>
#include <distributions/rintnunif.hpp>
#include <distributions/tstudent.hpp>
#include <distributions/wishart.hpp>

// [[Rcpp::depends("RcppArmadillo")]]

/* 
  Clean parameter, discard the middle not used values for the clusters and 
    update the correspondent parameters. 

  args:
    - Lambda:  array, each slice is a precision matrix
    - mu:      matrix, each row is a mean vector
    - clust:   vector, each (integer) value is the cluster of the corresp. obs.
    - useful:  vector, binary. Value 1 implies that the corresponding position is
                              used, value 0 implies that is not used

  Void function.
*/

void para_cleanser(arma::cube &Lambda, 
                   arma::mat &mu, 
                   arma::vec &clust, 
                   arma::vec &useful) {
  int k = (int) arma::sum(useful);
  
  // for all the used parameters
  for(int i = 0; i < k; i++){
    
    // if a cluster is empty
    if((int) arma::sum(clust == i) == 0){
      
      // find the last full cluster, then swap
      for(int j = k; j > i; j--){
        if((int) arma::sum(clust == j) != 0){
          clust( arma::find(clust == j) ).fill(i);
          mu.swap_rows(i,j);
          Lambda.slice(i).swap(Lambda.slice(j));
          break;
        }
      }
    }
    if(arma::sum(clust == i) == 0) useful[i] = 0;
  }
}
