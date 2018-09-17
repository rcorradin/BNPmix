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
#include "distributions.hpp"
// [[Rcpp::depends("RcppArmadillo")]]

/*
Frequencies

agrs:
- vec
return:
- freq
*/

//[[Rcpp::export]]
arma::vec freq_vec(arma::vec vector){
  arma::vec uniq = unique(vector);
  int n_item = uniq.n_elem;
  arma::vec result(n_item);
  for(int j = 0; j < n_item; j++){
    result[j] = (int) arma::sum(vector == uniq[j]);
  }
  return(result);
}

/*
 *  Eval density
 *  
 */

//[[Rcpp::export]]
arma::vec eval_density(arma::vec grid,
                       arma::vec mu,
                       arma::vec s2,
                       arma::vec probs){
  probs = probs / sum(probs);
  arma::vec result(grid.n_elem, arma::fill::zeros);
  for(int j = 0; j < mu.n_elem; j++){
    result += probs[j] * normpdf(grid, mu[j], sqrt(s2[j]));
  }
  return(result);
}

/*
 *  Eval density - MV
 *  
 */

//[[Rcpp::export]]
arma::vec eval_density_mv(arma::mat grid,
                          arma::mat mu,
                          arma::cube s2,
                          arma::vec probs){
  double d = (double) grid.n_cols;
  probs = probs / sum(probs);
  arma::vec result(grid.n_rows, arma::fill::zeros);
  for(int j = 0; j < probs.n_elem; j++){
    arma::mat rooti  = arma::trans(arma::inv(trimatu(arma::chol(s2.slice(j)))));
    double sum_rooti = arma::sum(log(rooti.diag()));
    for(int i = 0; i < grid.n_rows; i++){
      arma::vec cdata  = rooti * arma::trans(grid.row(i) - mu.row(j)) ;
      double out       = - (d / 2.0) * log(2.0 * M_PI) - 0.5 * arma::sum(cdata%cdata) + sum_rooti;
      result[i] += probs[j] * exp(out);
    }
  }
  return(result);
}
