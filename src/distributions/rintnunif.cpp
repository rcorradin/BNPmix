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
// [[Rcpp::depends("RcppArmadillo")]]

/* 
  generate and integer value from a
  non-uniform distribution
  args:
  	- p: vector of probabilities
  	- a: number of different values
  return an integer value
*/

int rintnunif(arma::vec prob, 
              int a){ 
  double u = arma::randu(1)[0];
  prob = arma::cumsum(prob);
  
  for(int k = 0; k < a; k++) {
    if(u <= prob[k]) {
      return k; 
    }
  }
  return 0;
}