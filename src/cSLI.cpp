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
#include "common_utilities.hpp"
#include "cSLI_functions.hpp"
// [[Rcpp::depends("RcppArmadillo")]]

/*
 *   cSLI - Univariate Conditional Slice sampler 
 *   
 *   args:
 *   - data:      a vector of observations
 *   - grid:      vector to evaluate the density
 *   - niter:     number of iterations
 *   - nburn:     number of burn-in iterations
 *   - m0:        expectation of location component
 *   - k0:        tuning parameter of variance of location component
 *   - a0, b0:    parameters of scale component
 *   - mass:      mass of Dirichlet process
 *   - nupd:      number of iterations to show current updating
 *                (default niter/10)
 *   - out_param: if TRUE, return also the location and scale paramteres lists
 *                (default FALSE)
 *   - out_dens:  if TRUE, return also the estimated density (default TRUE)
 *   - sigma_PY:   second parameter PY
 *    
 *   output, list:
 *   - dens:  matrix, each row a density evaluated on the grid
 *   - clust: matrix, each row an allocation
 *   - mu:    list, each element a locations vector
 *   - s2:    list, each element a scale vector
 *   - probs: list, each element a probabilities vector
 */
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
                int process = 0,
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
  
  int start_s = clock();
  int current_s;
  // strarting loop
  if(process== 0){
    for(int iter = 0; iter < niter; iter++){
      
      para_clean_SLI(mu, 
                     s2, 
                     clust,
                     v,
                     w);
      
      accelerate_SLI(data,
                     mu,
                     s2,
                     v,
                     w,
                     clust,
                     m0,
                     k0,
                     a0,
                     b0,
                     mass);
      
      int old_length = mu.n_elem;
      
      update_u_SLI(clust,
                   w,
                   u);
      
      grow_param_SLI(mu,
                     s2,
                     v,
                     w,
                     u,
                     m0,
                     k0,
                     a0,
                     b0,
                     mass,
                     n);
      
      update_cluster_SLI(data,
                         mu,
                         s2,
                         clust,
                         w,
                         u,
                         old_length,
                         iter,
                         new_val);
      
      
      if(iter >= nburn){
        
        result_clust.row(iter - nburn) = arma::trans(clust);
        result_mu.push_back(mu);
        result_s2.push_back(s2);
        result_probs.push_back(w);
        n_clust(iter - nburn) = w.n_elem;
        if(out_dens){
          result_dens.row(iter - nburn) = arma::trans(eval_density(grid,
                                                      mu,
                                                      s2,
                                                      w));
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
  } else if(process == 1){
    for(int iter = 0; iter < niter; iter++){
      
      para_clean_SLI(mu, 
                     s2, 
                     clust,
                     v,
                     w);
      
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
      
      update_u_SLI(clust,
                   w,
                   u);
      
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
      
      update_cluster_SLI(data,
                         mu,
                         s2,
                         clust,
                         w,
                         u,
                         old_length,
                         iter,
                         new_val);
      
      
      if(iter >= nburn){
        
        result_clust.row(iter - nburn) = arma::trans(clust);
        result_mu.push_back(mu);
        result_s2.push_back(s2);
        result_probs.push_back(w);
        n_clust(iter - nburn) = w.n_elem;
        if(out_dens){
          result_dens.row(iter - nburn) = arma::trans(eval_density(grid,
                                                      mu,
                                                      s2,
                                                      w));
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
    Rcpp::Rcout << "\n" << "WoooW! Have you seen how fast am I?\n";
  }
  
  Rcpp::List resu;
  if(out_param){
    resu["dens"]   = result_dens;
    resu["clust"]  = result_clust;
    resu["mu"]     = result_mu;
    resu["s2"]     = result_s2;
    resu["probs"]  = result_probs;
    resu["newval"] = new_val;
    resu["nclust"] = n_clust;
    resu["time"]   = double(end_s-start_s)/CLOCKS_PER_SEC;
  } else {
    resu["dens"]   = result_dens;
    resu["clust"]  = result_clust;
    resu["newval"] = new_val;
    resu["nclust"] = n_clust;
    resu["time"]   = double(end_s-start_s)/CLOCKS_PER_SEC;
  }
  return resu;
}


/*
*   cSLI_mv - MULTIVARIATE Conditional Slice sampler 
*   
*   args:
*   - data:      a vector of observations
*   - grid:      vector to evaluate the density
*   - niter:     number of iterations
*   - nburn:     number of burn-in iterations
*   - m0:        expectation of location component
*   - k0:        tuning parameter of variance of location component
*   - a0, b0:    parameters of scale component
*   - mass:      mass of Dirichlet process
*   - nupd:      number of iterations to show current updating
*                (default niter/10)
*   - out_param: if TRUE, return also the location and scale paramteres lists
*                (default FALSE)
*   - out_dens:  if TRUE, return also the estimated density (default TRUE)
*   - process:   if 0 DP, if 1 PY
*   - sigma_PY:  second parameter of PY   
*    
*   output, list:
*   - dens:  matrix, each row a density evaluated on the grid
*   - clust: matrix, each row an allocation
*   - mu:    list, each element a locations vector
*   - s2:    list, each element a scale vector
*   - probs: list, each element a probabilities vectorarma::vec new_val(niter);
*/
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
                   bool print_message = 1){
  
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
  arma::mat result_dens(niter - nburn, grid.n_rows);
  
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

  int start_s = clock();
  int current_s;
  // strarting loop
  if(process == 0){
    for(int iter = 0; iter < niter; iter++){
      
      para_clean_SLI_mv(mu, 
                        s2, 
                        clust,
                        v,
                        w);
      
      accelerate_SLI_mv(data,
                        mu,
                        s2,
                        v,
                        w,
                        clust,
                        m0,
                        k0,
                        S0,
                        n0,
                        mass);
      
      int old_length = mu.n_elem;
      
      update_u_SLI(clust,
                   w,
                   u);
      
      grow_param_SLI_mv(mu,
                        s2,
                        v,
                        w,
                        u,
                        m0,
                        k0,
                        S0,
                        n0,
                        mass,
                        n);
      
      update_cluster_SLI_mv(data,
                            mu,
                            s2,
                            clust,
                            w,
                            u,
                            old_length,
                            iter,
                            new_val);
      
      
      if(iter >= nburn){
        
        result_clust.row(iter - nburn) = arma::trans(clust);
        result_mu.push_back(mu);
        result_s2.push_back(s2);
        result_probs.push_back(w);
        n_clust(iter - nburn) = w.n_elem;
        if(out_dens){
          dens = eval_density_mv(grid,
                                 mu,
                                 s2,
                                 w);
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
  } else if(process == 1){
    for(int iter = 0; iter < niter; iter++){
      
      para_clean_SLI_mv(mu, 
                        s2, 
                        clust,
                        v,
                        w);
      
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
      
      update_u_SLI(clust,
                   w,
                   u);
      
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
      
      update_cluster_SLI_mv(data,
                            mu,
                            s2,
                            clust,
                            w,
                            u,
                            old_length,
                            iter,
                            new_val);
      
      
      if(iter >= nburn){
        
        result_clust.row(iter - nburn) = arma::trans(clust);
        result_mu.push_back(mu);
        result_s2.push_back(s2);
        result_probs.push_back(w);
        n_clust(iter - nburn) = w.n_elem;
        if(out_dens){
          dens = eval_density_mv(grid,
                                 mu,
                                 s2,
                                 w);
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
    Rcpp::Rcout << "\n" << "WoooW! Have you seen how fast am I?\n";
  }
  
  Rcpp::List resu;
  if(out_param){
    resu["dens"]   = result_dens;
    resu["clust"]  = result_clust;
    resu["mu"]     = result_mu;
    resu["s2"]     = result_s2;
    resu["probs"]  = result_probs;
    resu["newval"] = new_val;
    resu["nclust"] = n_clust;
    resu["time"]   = double(end_s-start_s)/CLOCKS_PER_SEC;
  } else {
    resu["dens"]   = result_dens;
    resu["clust"]  = result_clust;
    resu["newval"] = new_val;
    resu["nclust"] = n_clust;
    resu["time"]   = double(end_s-start_s)/CLOCKS_PER_SEC;
  }
  return resu;
}