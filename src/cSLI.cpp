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
#include "Distributions.h"
#include "CommonUtilities.h"
#include "SliFunctions.h"
// [[Rcpp::depends("RcppArmadillo")]]

//' @export cSLI
//' @name cSLI
//' @title C++ function to estimate Pitman-Yor univariate mixtures via slice sampler
//'
//'
//' @param data a vector of observations
//' @param grid vector to evaluate the density
//' @param niter number of iterations
//' @param nburn number of burn-in iterations
//' @param m0 expectation of location component
//' @param k0 tuning parameter of variance of location component
//' @param a0 parameter of scale component
//' @param b0 parameter of scale component
//' @param mass mass of Dirichlet process
//' @param nupd number of iterations to show current updating
//' @param out_param if TRUE, return also the location and scale paramteres lists
//' @param out_dens if TRUE, return also the estimated density (default TRUE)
//' @param sigma_PY second parameter of PY
//' @param print_message print the status

//[[Rcpp::export]]
Rcpp::List cSLI(arma::vec data,
                arma::vec grid,
                int niter,
                int nburn,
                double m0,
                double k0,
                double a0,
                double b0,
                double mass,
                int nupd = 0,
                bool out_param = 0,
                bool out_dens = 1,
                double sigma_PY = 0,
                bool print_message = 1){

  if(nupd == 0){
    nupd = int(niter / 10);
  }

  int n = data.n_elem;

  arma::mat result_clust(niter - nburn, n);
  std::list<arma::vec> result_mu;
  std::list<arma::vec> result_s2;
  std::list<arma::vec> result_probs;
  arma::mat result_dens(niter - nburn, grid.n_elem);
  arma::vec new_val(niter);
  arma::vec n_clust(niter - nburn);
  result_dens.fill(0);

  arma::vec clust(n);
  clust.fill(0);

  arma::vec mu(1);
  arma::vec s2(1);
  arma::vec v(1);
  arma::vec w(1);
  arma::vec u(n, arma::fill::randu);

  new_val.fill(0);
  n_clust.fill(0);
  mu.fill(m0);
  s2.fill(b0 / (a0 - 1));

  // time quantities
  int start_s = clock();
  int current_s;

  // strarting loop
  for(arma::uword iter = 0; iter < niter; iter++){

    // update the parameters
    accelerate_SLI_PY(data,
                      mu,
                      s2,
                      v,
                      w,
                      clust,
                      m0,
                      k0,
                      a0,
                      b0,
                      mass,
                      sigma_PY);

    int old_length = mu.n_elem;

    // update the stick breaking weights
    update_u_SLI(clust,
                 w,
                 u);

    // extend the stick breaking representation
    grow_param_SLI_PY(mu,
                      s2,
                      v,
                      w,
                      u,
                      m0,
                      k0,
                      a0,
                      b0,
                      mass,
                      n,
                      sigma_PY);

    // update the allocation
    update_cluster_SLI(data,
                       mu,
                       s2,
                       clust,
                       w,
                       u,
                       old_length,
                       iter,
                       new_val);

    // save the results
    if(iter >= nburn){
      result_mu.push_back(mu);
      result_s2.push_back(s2);
      result_probs.push_back(w);
      n_clust(iter - nburn) = w.n_elem;
      if(out_dens){
        result_dens.row(iter - nburn) = arma::trans(eval_density(grid, mu, s2, w));
      }
    }

    // clean the parameters
    para_clean_SLI(mu,
                   s2,
                   clust,
                   v,
                   w);

    // save the results
    if(iter >= nburn){
      result_clust.row(iter - nburn) = arma::trans(clust);
    }

    if(print_message){
      // print the current completed work
      if((iter + 1) % nupd == 0){
        current_s = clock();
        Rcpp::Rcout << "Completed:\t" << (iter + 1) << "/" << niter << " - in " <<
          double(current_s-start_s)/CLOCKS_PER_SEC << " sec\n";
      }
    }
    Rcpp::checkUserInterrupt();
  }

  int end_s = clock();
  if(print_message){
    Rcpp::Rcout << "\n" << "Estimation done in " << double(end_s-start_s)/CLOCKS_PER_SEC << " seconds\n";
  }

  Rcpp::List resu;
  if(out_param){
    resu["dens"]   = result_dens;
    resu["clust"]  = result_clust;
    resu["mu"]     = result_mu;
    resu["s2"]     = result_s2;
    resu["probs"]  = result_probs;
    resu["newval"] = new_val;
    resu["time"]   = double(end_s-start_s)/CLOCKS_PER_SEC;
  } else {
    resu["dens"]   = result_dens;
    resu["clust"]  = result_clust;
    resu["newval"] = new_val;
    resu["time"]   = double(end_s-start_s)/CLOCKS_PER_SEC;
  }
  return resu;
}

//' @export cSLI_mv
//' @name cSLI_mv
//' @title C++ function to estimate Pitman-Yor multivariate mixtures via slice sampler
//'
//'
//' @param data a matrix of observations
//' @param grid matrix of points to evaluate the density
//' @param niter number of iterations
//' @param nburn number of burn-in iterations
//' @param m0 expectation of location component
//' @param k0 tuning parameter of variance of location component
//' @param S0 parameter of scale component
//' @param n0 parameter of scale component
//' @param mass mass of Dirichlet process
//' @param nupd number of iterations to show current updating
//' @param out_param if TRUE, return also the location and scale paramteres lists
//' @param out_dens if TRUE, return also the estimated density (default TRUE)
//' @param sigma_PY second parameter of PY
//' @param print_message print the status
//' @param light_dens if TRUE return only the posterior mean of the density

//[[Rcpp::export]]
Rcpp::List cSLI_mv(arma::mat data,
                   arma::mat grid,
                   int niter,
                   int nburn,
                   arma::vec m0,
                   double k0,
                   arma::mat S0,
                   double n0,
                   double mass,
                   int nupd = 0,
                   bool out_param = 0,
                   bool out_dens = 1,
                   int process = 0,
                   double sigma_PY = 0,
                   bool print_message = 1,
                   bool light_dens = 1){

  if(nupd == 0){
    nupd = int(niter / 10);
  }

  int n = data.n_rows;
  int d = data.n_cols;

  // initialize results objects
  arma::mat result_clust(niter - nburn, n);
  std::list<arma::mat> result_mu;
  std::list<arma::cube> result_s2;
  std::list<arma::vec> result_probs;
  arma::mat result_dens(grid.n_rows, 1);
  if(!light_dens){
    result_dens.resize(niter - nburn, grid.n_rows);
  }
  result_dens.fill(0);

  // initialize required object inside the loop
  arma::vec clust(n);
  arma::mat mu(1,d);
  arma::cube s2(d,d,1);
  arma::vec v(1);
  arma::vec w(1);
  arma::vec u(n, arma::fill::randu);
  arma::vec dens(grid.n_rows);
  arma::vec new_val(niter);
  arma::vec n_clust(niter - nburn);

  // fill the initialized quantity
  clust.fill(0);
  mu.row(0) = arma::trans(m0);
  s2.slice(0) = S0 / (n0 - d - 1);
  new_val.fill(0);
  n_clust.fill(0);

  // time quantities
  int start_s = clock();
  int current_s;

  // strarting loop
  for(arma::uword iter = 0; iter < niter; iter++){

    // update the parameters
    accelerate_SLI_PY_mv(data,
                         mu,
                         s2,
                         v,
                         w,
                         clust,
                         m0,
                         k0,
                         S0,
                         n0,
                         mass,
                         sigma_PY);

    int old_length = mu.n_elem;

    // update the stick breaking weights
    update_u_SLI(clust,
                 w,
                 u);

    // extend the stick breaking representation
    grow_param_SLI_PY_mv(mu,
                         s2,
                         v,
                         w,
                         u,
                         m0,
                         k0,
                         S0,
                         n0,
                         mass,
                         n,
                         sigma_PY);

    // update the allocation
    update_cluster_SLI_mv(data,
                          mu,
                          s2,
                          clust,
                          w,
                          u,
                          old_length,
                          iter,
                          new_val);

    // save the results
    if(iter >= nburn){
      result_mu.push_back(mu);
      result_s2.push_back(s2);
      result_probs.push_back(w);
      n_clust(iter - nburn) = w.n_elem;
      if(out_dens){
        dens = eval_density_mv(grid,
                               mu,
                               s2,
                               w);
        if(light_dens){
          result_dens += dens;
        } else {
          result_dens.row(iter - nburn) = arma::trans(dens);
        }
      }
    }

    // clean the parameters
    para_clean_SLI_mv(mu,
                      s2,
                      clust,
                      v,
                      w);

    // save the results
    if(iter >= nburn){
      result_clust.row(iter - nburn) = arma::trans(clust);
    }

    if(print_message){
      // print the current completed work
      if((iter + 1) % nupd == 0){
        current_s = clock();
        Rcpp::Rcout << "Completed:\t" << (iter + 1) << "/" << niter << " - in " <<
          double(current_s-start_s)/CLOCKS_PER_SEC << " sec\n";
      }
    }
    Rcpp::checkUserInterrupt();
  }

  int end_s = clock();
  if(print_message){
    Rcpp::Rcout << "\n" << "Estimation done in " << double(end_s-start_s)/CLOCKS_PER_SEC << " seconds\n";
  }

  Rcpp::List resu;
  if(out_param){
    if(light_dens){
      resu["dens"]   = result_dens / (niter - nburn);
    } else {
      resu["dens"]   = result_dens;
    }
    resu["clust"]  = result_clust;
    resu["mu"]     = result_mu;
    resu["s2"]     = result_s2;
    resu["probs"]  = result_probs;
    resu["newval"] = new_val;
    resu["time"]   = double(end_s-start_s)/CLOCKS_PER_SEC;
  } else {
    if(light_dens){
      resu["dens"]   = result_dens / (niter - nburn);
    } else {
      resu["dens"]   = result_dens;
    }
    resu["clust"]  = result_clust;
    resu["newval"] = new_val;
    resu["time"]   = double(end_s-start_s)/CLOCKS_PER_SEC;
  }
  return resu;
}
