/*==================================================================================
 Copyright (C) 2019 Riccardo Corradin

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
#include <cmath>
// [[Rcpp::depends("RcppArmadillo")]]

//' @export BNPmix_psm
//' @name BNPmix_psm
//' @title C++ function - compute the posterior similarity matrix
//' @keywords internal
//'
//' @param M a matrix (r x n), r number of replications, n number of observations
//'
//' @examples{
//'   M <- matrix(c(1,1,1,2,1,1,2,2,1,1,2,1,1,2,1,1), ncol = 4)
//'   BNPmix_psm(M)
//' }
//'

//[[Rcpp::export]]
arma::mat BNPmix_psm(arma::mat M){
  // initialize results
  arma::mat result(M.n_cols, M.n_cols, arma::fill::zeros);

  for(arma::uword i = 0; i < M.n_cols; i++){
    for(arma::uword j = 0; j <= i; j++){
      result(i,j) = arma::accu(M.col(i) == M.col(j));
      result(j,i) = result(i,j);
    }
    Rcpp::checkUserInterrupt();
  }
  return(result / M.n_rows);
}

//' @export clean_partition
//' @name clean_partition
//' @title C++ function - clean the partition matrix
//' @keywords internal
//'
//' @param M a matrix (r x n), r number of replications, n number of observations
//'
//' @examples{
//'   M <- matrix(c(1,1,1,3,1,1,4,4,1,1,3,1,1,3,1,1), ncol = 4)
//'   clean_partition(M)
//' }
//'

//[[Rcpp::export]]
arma::mat clean_partition(arma::mat M){

  arma::uvec index(M.n_cols);
  arma::vec tvec(M.n_cols);
  // initialize results
  arma::mat result(M.n_rows, M.n_cols, arma::fill::zeros);

  // for each row
  for(arma::uword k = 0; k < M.n_rows; k++){
    tvec = M.row(k).t();

    for(arma::uword j = 0; j < max(M.row(k)); j++){
      while((arma::accu(tvec == j + 1) == 0) && (arma::accu(tvec > j + 1) != 0)){
        index = find(tvec > j + 1);
        tvec(index) = tvec(index) - 1;
      }
    }

    result.row(k) = tvec.t();
    Rcpp::checkUserInterrupt();
  }
  return(result);
}

//' @export BNPmix_VI_LB
//' @name BNPmix_VI_LB
//' @title C++ function - compute the VI lower bound
//' @keywords internal
//'
//' @param M a matrix (r x n), r number of replications, n number of observations
//' @param psm_mat a posterior similarity matrix
//'
//' @examples{
//'   M <- matrix(c(1,1,1,2,1,1,2,2,1,1,2,1,1,1,1,2), ncol = 4)
//'   psmM <- BNPmix_psm(M)
//'   BNPmix_VI_LB(M, psmM)
//' }
//'

//[[Rcpp::export]]
arma::vec BNPmix_VI_LB(arma::mat C_mat, arma::mat psm_mat){

  arma::vec result(C_mat.n_rows);
  double f = 0.0;
  int n = psm_mat.n_cols;
  arma::vec tvec(n);

  for(arma::uword j = 0; j < C_mat.n_rows; j++){
    f = 0.0;
    for(arma::uword i = 0; i < n; i++){
      tvec = psm_mat.col(i);
      f += (log2(arma::accu(C_mat.row(j) == C_mat(j,i))) +
        log2(arma::accu(tvec)) -
        2 * log2(arma::accu(tvec.elem(arma::find(C_mat.row(j).t() == C_mat(j,i))))))/n;
    }
    result(j) = f;
    Rcpp::checkUserInterrupt();
  }
  return(result);
}

//' @export BNPmix_BIN
//' @name BNPmix_BIN
//' @title C++ function - compute the Binder distances
//' @keywords internal
//'
//' @param M a matrix (r x n), r number of replications, n number of observations
//' @param psm_mat a posterior similarity matrix
//'
//' @examples{
//'   M <- matrix(c(1,1,1,2,1,1,2,2,1,1,2,1,1,1,1,2), ncol = 4)
//'   psmM <- BNPmix_psm(M)
//'   BNPmix_BIN(M, psmM)
//' }
//'

//[[Rcpp::export]]
arma::vec BNPmix_BIN(arma::mat C_mat, arma::mat psm_mat){

  arma::vec result(C_mat.n_rows);
  arma::mat tmat(psm_mat);
  int n = C_mat.n_cols;

  for(arma::uword j = 0; j < C_mat.n_rows; j++){

    tmat.fill(0.0);
    for(arma::uword i = 0; i < n; i++){
      for(arma::uword k = 0; k < n; k++){
        if(C_mat(j,i) == C_mat(j,k)){ tmat(i,k) = 1; }
      }
    }
    result(j) = arma::accu(abs(tmat - psm_mat))/2;
    Rcpp::checkUserInterrupt();
  }
  return(result);
}
