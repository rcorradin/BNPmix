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
#include "DDP_functions.hpp"
// [[Rcpp::depends("RcppArmadillo")]]

/*
*   cPUS - UNIVARIATE Conditional dependent Dirichlet process sampler (DDP)
*   
*   args:
*   - data:      a vector of observations
*   - group:     a vector of groups
*   - grid:      vector to evaluate the density
*   - niter:     number of iterations
*   - nburn:     number of burn-in iterations
*   - m0:        expectation of location component
*   - k0:        tuning parameter of variance of location component
*   - a0, b0:    parameters of scale component
*   - mass:      mass of Dirichlet process
*   - napprox:   number of approximating values
*   - nupd:      number of iterations to show current updating
*                (default niter/10)
*   - out_dens:  if TRUE, return also the estimated density (default TRUE)
*   - process:   if 0 DP, if 1 PY
*   - sigma_PY:  second parameter of PY
*    
*   output, list:
*   - dens:  matrix, each row a density evaluated on the grid
*   - clust: matrix, each row an allocation
*   - mu:    list, each element a locations vector
*   - s2:    list, each element a scale vector
*   - probs: list, each element a probabilities vector
*/

//[[Rcpp::export]]
Rcpp::List cDDP(arma::vec data,
                arma::vec group,
                int ngr,
                arma::vec grid,
                int niter,
                int nburn,
                double m0,
                double k0,
                double a0,
                double b0,
                double mass,
                double wei,
                int napprox,
                int n_approx_unif,
                int nupd = 0,
                bool out_dens = 1,
                bool print_message = 1){
  
  // 
  // INITIALIZE
  // 
  
  if(nupd == 0){
    nupd = (int) (niter / 10);
  }
  int n = data.n_elem;
  
  // initialize results objects 
  arma::mat result_clust(niter - nburn, n);
  arma::mat result_zeta(niter - nburn, n);
  arma::cube result_dens(grid.n_elem, ngr, niter - nburn);
  
  // initialize required object inside the loop
  arma::vec clust(n);
  arma::vec zeta(group);
  arma::field<arma::vec> mu(ngr + 1);
  arma::field<arma::vec> s2(ngr + 1);
  arma::field<arma::vec> mutemp(ngr + 1);
  arma::field<arma::vec> s2temp(ngr + 1);
  arma::field<arma::vec> freqtemp(ngr + 1);
  arma::field<arma::vec> ptilde(ngr + 1);
  arma::mat dens(grid.n_elem, ngr);
  arma::vec new_val(niter);
  arma::mat n_clust(niter - nburn, ngr + 1);
  arma::field<arma::vec> mujoin(ngr + 1);
  arma::field<arma::vec> s2join(ngr + 1);
  arma::field<arma::vec> probjoin(ngr + 1);
  arma::vec max_val(ngr + 1);
  arma::vec w(ngr);
  
  // fill the initialized quantity
  for(arma::uword j = 0; j <= ngr; j++){
    
    // resize
    mu(j).resize(1);
    s2(j).resize(1);
    mutemp(j).resize(1);
    s2temp(j).resize(1);
    freqtemp(j).resize(1);
    ptilde(j).resize(2);
    
    // fill
    mu(j).fill(m0);
    s2(j).fill(b0 / (a0 - 1));
    mutemp(j).fill(0);
    s2temp(j).fill(0);
    freqtemp(j).fill(1);
    ptilde(j) = rdirich_mass(freq_vec(zeta(zeta == j)), mass * wei);
  }
  
  ptilde(0).resize(1);
  ptilde(0).fill(1);
  
  new_val.fill(0);
  n_clust.fill(0);
  clust.fill(0);
  
  
  int start_s = clock();
  int current_s;
  
  //
  // STARTING
  // MAIN 
  // LOOP
  // 
  
  for(int iter = 0; iter < niter; iter++){
    
    // clean parameter objects
    para_clean_DDP(mu,
                   s2,
                   clust,
                   group,
                   zeta,
                   ngr);
    
    // sample the probability vector
    // from Dirichlet distribution
    for(arma::uword g = 0; g <= ngr; g++){
      if(any(zeta == g)){
        if(g != 0){
          ptilde(g) = rdirich_mass(freq_vec(clust.elem(arma::find(zeta == g))), mass * wei);    
        } else {
          ptilde(g) = rdirich_mass(freq_vec(clust.elem(arma::find(zeta == g))), mass * (1 - wei));  
        }
        
      } else {
        ptilde(g).resize(1);
        ptilde(g).fill(1);
      }
      
    }
    
    // update w
    update_w_DDP(w,
                 mass,
                 wei,
                 n_approx_unif,
                 zeta,
                 clust,
                 group,
                 ngr);
    
    // update zeta
    zeta_update_DDP(zeta,
                    group,
                    clust,
                    mu,
                    s2,
                    ptilde,
                    w,
                    mass,
                    wei,
                    ngr);
    
    // clean parameter objects
    para_clean_DDP(mu,
                   s2,
                   clust,
                   group,
                   zeta,
                   ngr);
    
    // acceleration step
    accelerate_DDP(data,
                   group,
                   zeta,
                   mu,
                   s2,
                   clust,
                   m0,
                   k0,
                   a0,
                   b0,
                   ngr);
    
    // sample the probability vector
    // from Dirichlet distribution
    for(arma::uword g = 0; g <= ngr; g++){
      if(any(zeta == g)){
        if(g != 0){
          ptilde(g) = rdirich_mass(freq_vec(clust.elem(arma::find(zeta == g))), mass * wei);    
        } else {
          ptilde(g) = rdirich_mass(freq_vec(clust.elem(arma::find(zeta == g))), mass * (1 - wei));  
        }
        
      } else {
        ptilde(g).resize(1);
        ptilde(g).fill(1);
      }
      
    }

    // simulate the required values
    simu_trunc_DDP(mutemp,
                   s2temp,
                   freqtemp,
                   mass,
                   wei,
                   m0,
                   k0,
                   a0,
                   b0,
                   napprox,
                   ngr);
    
    // joint the existent parameters with
    // the simulated new ones
    for(arma::uword g = 0; g <= ngr; g++){
      
      if(any(zeta == g)){
        ptilde(g) = rdirich_mass(freq_vec(clust.elem(arma::find(zeta == g))), mass);  
        
        mujoin(g) = arma::join_cols(mu(g), mutemp(g));
        s2join(g) = arma::join_cols(s2(g), s2temp(g));
        int nkpt = ptilde(g).n_elem;
        arma::vec index = arma::regspace(0, nkpt - 1);
        probjoin(g) = arma::join_cols(ptilde(g).elem(arma::find(index < nkpt - 1)),
                 ptilde(g)[nkpt - 1] * freqtemp(g) / napprox);
        max_val(g) = mu(g).n_elem;  
        
      } else {
        ptilde(g) = 1.0;
        
        mujoin(g) = mutemp(g);
        s2join(g) = s2temp(g);
        probjoin(g) = freqtemp(g) / napprox;
        max_val(g)  = 0;  
        
      }
    }
    
    // update cluster allocations
    clust_update_DDP(data,
                     group,
                     zeta,
                     mujoin,
                     s2join,
                     probjoin,
                     clust,
                     max_val,
                     iter,
                     ngr,
                     new_val);
    
    mu = mujoin;
    s2 = s2join;
    
    // if the burn-in phase is complete
    if(iter >= nburn){
      
      // save clust and zeta in the export objects
      result_clust.row(iter - nburn) = arma::trans(arma::conv_to<arma::vec>::from(clust));
      result_zeta.row(iter - nburn)  = arma::trans(arma::conv_to<arma::vec>::from(zeta));
      
      // compute the density for each urn
      for(arma::uword g = 1; g <= ngr; g++){
        // n_clust(iter - nburn, g) = probjoin.n_elem;
        if(out_dens){
          
          dens.col(g-1) = w(g - 1) * eval_density(grid, mujoin(g),
                   s2join(g), probjoin(g)) + (1 - w(g - 1)) * eval_density(grid, mujoin(0),
                   s2join(0), probjoin(0));
          
          // dens.col(g-1) = sum(group == g && zeta == g) / sum(group == g) * eval_density(grid, mujoin(g),
          //          s2join(g), probjoin(g)) + (1 - sum(group == g && zeta == g) / sum(group == g)) * 
          //            eval_density(grid, mujoin(0), s2join(0), probjoin(0));
        }
      }
      
      // update the slice of the current iteration in the density object
      result_dens.slice(iter - nburn) = dens;
    }

    if(print_message){
      // print the current completed work
      if((iter + 1) % nupd == 0){
        current_s = clock();
        Rcpp::Rcout << "Completed:\t" << (iter + 1) << "/" << niter << " - in " <<
          double(current_s-start_s)/CLOCKS_PER_SEC << " sec\n";
      }
    }
    
    // check for aborting comand
    Rcpp::checkUserInterrupt();
  }
  int end_s = clock();
  if(print_message){
    Rcpp::Rcout << "\n" << "WoooW! Have you seen how fast am I?\n";
  }
  
  // export a list with the results
  Rcpp::List resu;
  resu["dens"]   = result_dens;
  resu["clust"]  = result_clust;
  resu["zeta"]   = result_zeta;
  resu["newval"] = new_val;
  resu["nclust"] = n_clust;
  resu["time"]   = double(end_s-start_s)/CLOCKS_PER_SEC;
  
  return resu;
}