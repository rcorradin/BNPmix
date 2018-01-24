#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

//----------------------
// best partition finder
//----------------------

// [[Rcpp::export]]
Rcpp::List find_part(arma::mat clust_mat) {
  int nrow = clust_mat.n_rows;
  int n    = clust_mat.n_cols;
  arma::cube results(n,n,nrow);
  
  for(int i = 0; i < nrow; i++){
    for(int j = 0; j < n; j++){
      for(int k = j; k < n; k++){
        results(j,k,i) = results(k,j,i) = (clust_mat(i,j) == clust_mat(i,k) ? 1 : 0);
      }
    }
  }
  
  arma::mat av_mat = sum(results, 2) / nrow;
  arma::vec dist_v(nrow);
  for(int i = 0; i < nrow; i++){
    dist_v(i) = arma::accu(abs(av_mat - results.slice(i)));
  }
  
  arma::mat best = results.slice(arma::index_min(dist_v));
  
  Rcpp::List resu;
  resu["average"] = av_mat;
  resu["best"]    = best;
  return resu;
  
}

//--------------
// pairwise mat
//--------------

// [[Rcpp::export]]
Rcpp::List pairwise_mat(arma::mat clust_mat) {
  int nrow = clust_mat.n_rows;
  int n    = clust_mat.n_cols;
  arma::mat results(n,n);
  
  for(int i = 0; i < nrow; i++){
    for(int j = 0; j < n; j++){
      for(int k = j; k < n; k++){
        results(j,k) += (clust_mat(i,j) == clust_mat(i,k) ? 1 : 0);
        results(k,j) += (clust_mat(i,j) == clust_mat(i,k) ? 1 : 0);
      }
    }
  }
  
  for(int i = 0; i < n; i++){
    results(i,i) /= 2;
  }
  Rcpp::List resu;
  resu["PW"] = results;
  return resu;
}


//--------------
// estimated ISE
//--------------

// [[Rcpp::export]]
double est_ISE_2D(arma::vec estimated, arma::vec teoric, arma::mat grid) {
  arma::vec dim1 = arma::unique(grid.col(0));
  arma::vec dim2 = arma::unique(grid.col(1));
  return arma::sum(pow(estimated - teoric, 2)) * (dim1(1) - dim1(0)) * (dim2(1) - dim2(0));
}