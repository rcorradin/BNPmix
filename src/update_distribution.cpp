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

#include "gaussian.hpp"
#include "rintnunif.hpp"
#include "tstudent.hpp"
#include "wishart.hpp"

// [[Rcpp::depends("RcppArmadillo")]]

/*
  Update distribution (approximate).

  args:
    - grid:    matrix, given points to evaluate the density
    - Lambda:  array, each slice is a precision matrix
    - mu:      matrix, each row is a mean vector
    - clust:   vector, each (integer) value is the cluster of the corresp. obs.
    - useful:  vector, binary. Value 1 implies that the corresponding position is
                              used, value 0 implies that is not used
    - theta    double, precision parameter of the Dirichlet process
    - n:       int, number of observation

  Void function.
*/

arma::vec update_distribution(arma::mat grid,
                              int grid_l,
                              arma::mat mu,
                              arma::cube Lambda,
                              arma::vec useful,
                              arma::vec clust,
                              double theta,
                              int n){
  int k = (int) arma::sum(useful);
  arma::vec temp_out(grid_l);
  temp_out.fill(0);

  // for each different component
  for(int j = 0; j < k; j++){

    // evaluated the density weigthed by the frequence of the component
    temp_out += arma::sum(clust == j) * dmvnrm_ar_mat(grid, mu.row(j), Lambda.slice(j));
  }

  return temp_out / n;
}
