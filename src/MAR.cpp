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
#include "MarFunctions.h"
// [[Rcpp::depends("RcppArmadillo")]]

//' @export MAR
//' @name MAR
//' @title C++ function to estimate Pitman-Yor univariate mixtures via marginal sampler
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
//' @param process if 0 DP, if 1 PY
//' @param sigma_PY second parameter of PY
//' @param print_message print the status

//[[Rcpp::export]]
Rcpp::List MAR(arma::vec data,
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
              int process = 0,
              double sigma_PY = 0,
              bool print_message = 1){

  if(nupd == 0){
    nupd = (int) (niter / 10);
  }
  int n = data.n_elem;

  // initialize results objects
  arma::mat result_clust(niter - nburn, n);
  std::list<arma::vec> result_mu(niter - nburn);
  std::list<arma::vec> result_s2(niter - nburn);
  std::list<arma::vec> result_probs(niter - nburn);
  arma::mat result_dens(niter - nburn, grid.n_elem);

  // initialize required object inside the loop
  arma::vec clust(n);
  arma::vec mu(1);
  arma::vec s2(1);
  arma::vec dens(grid);
  arma::vec new_val(niter);
  arma::vec n_clust(niter - nburn);

  // fill the initialized quantity
  clust.fill(0);
  mu.fill(m0);
  s2.fill(b0 / (a0 - 1));
  new_val.fill(0);
  n_clust.fill(0);


  int start_s = clock();
  int current_s;
  // strarting loop
  if(process == 0){
    for(arma::uword iter = 0; iter < niter; iter++){

      // accelerating!!!
      accelerate_MAR(data,
                     mu,
                     s2,
                     clust,
                     m0,
                     k0,
                     a0,
                     b0);

      // update cluster allocation
      clust_update_MAR(data,
                       mu,
                       s2,
                       clust,
                       mass,
                       m0,
                       k0,
                       a0,
                       b0,
                       iter,
                       new_val);

      // clean parameter objects
      para_clean_MAR(mu,
                     s2,
                     clust);

      // if the burn-in phase is complete
      if(iter >= nburn){
        result_clust.row(iter - nburn) = arma::trans(arma::conv_to<arma::vec>::from(clust));
        result_mu.push_back(mu);
        result_s2.push_back(s2);
        arma::vec tab_freq = freq_vec(clust);
        result_probs.push_back(tab_freq);
        n_clust(iter - nburn) = tab_freq.n_elem;
        if(out_dens){
          dens = eval_density(grid,
                              mu,
                              s2,
                              tab_freq);
          result_dens.row(iter - nburn) = arma::trans(dens);
        }
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
  }
  else if(process == 1){
    for(arma::uword iter = 0; iter < niter; iter++){

      // accelerating!!!
      accelerate_MAR(data,
                     mu,
                     s2,
                     clust,
                     m0,
                     k0,
                     a0,
                     b0);

      // update cluster allocation
      clust_update_MAR_PY(data,
                          mu,
                          s2,
                          clust,
                          mass,
                          m0,
                          k0,
                          a0,
                          b0,
                          iter,
                          new_val,
                          sigma_PY);

      // clean parameter objects
      para_clean_MAR(mu,
                     s2,
                     clust);

      // if the burn-in phase is complete
      if(iter >= nburn){
        result_clust.row(iter - nburn) = arma::trans(arma::conv_to<arma::vec>::from(clust));
        result_mu.push_back(mu);
        result_s2.push_back(s2);
        arma::vec tab_freq = freq_vec(clust);
        result_probs.push_back(tab_freq);
        n_clust(iter - nburn) = tab_freq.n_elem;
        if(out_dens){
          dens = eval_density(grid,
                              mu,
                              s2,
                              tab_freq);
          result_dens.row(iter - nburn) = arma::trans(dens);
        }
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

//' @export MAR_mv
//' @name MAR_mv
//' @title C++ function to estimate Pitman-Yor multivariate mixtures via marginal sampler
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
//' @param process if 0 DP, if 1 PY
//' @param sigma_PY second parameter of PY
//' @param print_message print the status
//' @param light_dens if TRUE return only the posterior mean of the density

//[[Rcpp::export]]
Rcpp::List MAR_mv(arma::mat data,
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
    nupd = (int) (niter / 10);
  }
  int n = data.n_rows;
  int d = data.n_cols;

  // initialize results objects
  arma::mat result_clust(niter - nburn, n);
  std::list<arma::mat> result_mu(niter - nburn);
  std::list<arma::cube> result_s2(niter - nburn);
  std::list<arma::vec> result_probs(niter - nburn);
  arma::mat result_dens(grid.n_rows, 1);
  if(!light_dens){
    result_dens.resize(niter - nburn, grid.n_rows);
  }
  arma::vec n_clust(niter - nburn);

  // initialize required object inside the loop
  arma::vec clust(n);
  arma::mat mu(1,d);
  arma::cube s2(d,d,1);
  arma::vec dens(grid.n_rows);
  arma::vec new_val(niter);

  // fill the initialized quantity
  clust.fill(0);
  mu.row(0) = arma::trans(m0);
  s2.slice(0) = S0 / (n0 - d - 1);
  new_val.fill(0);
  n_clust.fill(0);

  int start_s = clock();
  int current_s;
  // strarting loop
  if(process == 0){
    for(arma::uword iter = 0; iter < niter; iter++){

      // accelerating
      accelerate_MAR_mv(data,
                        mu,
                        s2,
                        clust,
                        m0,
                        k0,
                        S0,
                        n0);

      // update cluster allocation
      clust_update_MAR_mv(data,
                          mu,
                          s2,
                          clust,
                          mass,
                          m0,
                          k0,
                          S0,
                          n0,
                          iter,
                          new_val);

      // clean parameter objects
      para_clean_MAR_mv(mu,
                        s2,
                        clust);

      // if the burn-in phase is complete
      if(iter >= nburn){
        result_clust.row(iter - nburn) = arma::trans(arma::conv_to<arma::vec>::from(clust));
        result_mu.push_back(mu);
        result_s2.push_back(s2);
        arma::vec tab_freq = freq_vec(clust);
        result_probs.push_back(tab_freq);
        n_clust(iter - nburn) = tab_freq.n_elem;
        if(out_dens){
          dens = eval_density_mv(grid,
                                 mu,
                                 s2,
                                 tab_freq);
          if(light_dens){
            result_dens += dens;
          } else {
            result_dens.row(iter - nburn) = arma::trans(dens);
          }
        }
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
  }
  else if(process == 1){
    for(arma::uword iter = 0; iter < niter; iter++){

      // accelerating
      accelerate_MAR_mv(data,
                        mu,
                        s2,
                        clust,
                        m0,
                        k0,
                        S0,
                        n0);

      // update cluster allocation
      clust_update_MAR_PY_mv(data,
                             mu,
                             s2,
                             clust,
                             mass,
                             m0,
                             k0,
                             S0,
                             n0,
                             iter,
                             new_val,
                             sigma_PY);

      // clean parameter objects
      para_clean_MAR_mv(mu,
                        s2,
                        clust);

      // if the burn-in phase is complete
      if(iter >= nburn){
        result_clust.row(iter - nburn) = arma::trans(arma::conv_to<arma::vec>::from(clust));
        result_mu.push_back(mu);
        result_s2.push_back(s2);
        arma::vec tab_freq = freq_vec(clust);
        result_probs.push_back(tab_freq);
        n_clust(iter - nburn) = tab_freq.n_elem;
        if(out_dens){
          dens = eval_density_mv(grid,
                                 mu,
                                 s2,
                                 tab_freq);
          if(light_dens){
            result_dens += dens;
          } else {
            result_dens.row(iter - nburn) = arma::trans(dens);
          }
        }
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
