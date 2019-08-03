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
#include "DdpFunctions.h"
// [[Rcpp::depends("RcppArmadillo")]]


//' @export
//' @name cDDP
//' @title C++ function to estimate DDP models with 1 grouping variables
//' @keywords internal
//'
//' @param data a vector of observations.
//' @param group group allocation of the data.
//' @param ngr number of groups.
//' @param grid vector to evaluate the density.
//' @param niter number of iterations.
//' @param nburn number of burn-in iterations.
//' @param m0 expectation of location component.
//' @param k0 tuning parameter of variance of location component.
//' @param a0 parameter of scale component.
//' @param b0 parameter of scale component.
//' @param mass mass of Dirichlet process.
//' @param wei prior weight of the specific processes.
//' @param b tuning parameter of weights distribution
//' @param napprox number of approximating values.
//' @param n_approx_unif number of approximating values of the importance step for the weights updating.
//' @param nupd number of iterations to show current updating.
//' @param out_dens if TRUE, return also the estimated density (default TRUE).
//' @param print_message print the status.
//' @param light_dens if TRUE return only the posterior mean of the density
//'


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
                bool print_message = 1,
                bool light_dens = 1){

  // INITIALIZE
  if(nupd == 0){
    nupd = (int) (niter / 10);
  }
  int n = data.n_elem;

  // initialize results objects
  arma::mat result_clust(niter - nburn, n);
  arma::mat result_group_log(niter - nburn, n);
  arma::cube result_dens(grid.n_elem, ngr, 1);
  if(!light_dens){
    result_dens.resize(grid.n_elem, ngr, niter - nburn);
  }
  arma::mat res_wei(niter - nburn, ngr);
  result_dens.fill(0);

  // initialize required object inside the loop
  arma::vec clust(n);
  arma::vec group_log(group);
  arma::field<arma::vec> mu(ngr + 1);
  arma::field<arma::vec> s2(ngr + 1);
  arma::field<arma::vec> mutemp(ngr + 1);
  arma::field<arma::vec> s2temp(ngr + 1);
  arma::field<arma::vec> freqtemp(ngr + 1);
  arma::field<arma::vec> ptilde(ngr + 1);
  arma::field<arma::vec> mujoin_complete(ngr);
  arma::field<arma::vec> s2join_complete(ngr);
  arma::field<arma::vec> probjoin_complete(ngr);
  arma::mat dens(grid.n_elem, ngr);
  arma::field<arma::vec> mujoin(ngr + 1);
  arma::field<arma::vec> s2join(ngr + 1);
  arma::field<arma::vec> probjoin(ngr + 1);
  arma::vec max_val(ngr);
  arma::vec w(ngr);
  arma::mat temp_proc_cum(n, 2);

  // fill the initialized quantity
  for(arma::uword j = 0; j <= ngr; j++){

    // resize
    mu(j).resize(1);
    s2(j).resize(1);
    mutemp(j).resize(1);
    s2temp(j).resize(1);
    freqtemp(j).resize(1);
    ptilde(j).resize(1);

    // fill
    mu(j).fill(m0);
    s2(j).fill(b0 / (a0 - 1));
    mutemp(j).fill(0);
    s2temp(j).fill(0);
    freqtemp(j).fill(1);
    ptilde(j).fill(1);
  }

  ptilde(0).resize(1);
  ptilde(0).fill(1);
  clust.fill(0);
  temp_proc_cum.fill(1.0);

  int start_s = clock();
  int current_s;

  //
  // STARTING
  // MAIN
  // LOOP
  //

  for(arma::uword iter = 0; iter < niter; iter++){

    //clean parameter objects
    para_clean_DDP(mu,
                    s2,
                    clust,
                    group,
                    group_log,
                    ngr);

    // sample the probability vector
    // from Dirichlet distribution
    for(arma::uword g = 0; g <= ngr; g++){
      if(any(group_log == g)){
        if(g != 0){
          ptilde(g) = rdirich_mass(freq_vec(clust.elem(arma::find(group_log == g))), mass * wei);
        } else {
          ptilde(g) = rdirich_mass(freq_vec(clust.elem(arma::find(group_log == g))), mass * (1 - wei));
        }
      } else {
        ptilde(g).resize(1);
        double temp = arma::randg(1, arma::distr_param(1.0, 1.0))[0];
        double temp2 = arma::randg(1, arma::distr_param(mass * wei, 1.0))[0];
        ptilde(g)(0) =  temp2 / ( temp + temp2);
      }

    }

    // update w
    update_w_DDP(w,
                  mass,
                  wei,
                  n_approx_unif,
                  group_log,
                  clust,
                  group,
                  temp_proc_cum,
                  ngr);

    // acceleration step
    accelerate_DDP(data,
                    group,
                    group_log,
                    mu,
                    s2,
                    clust,
                    m0,
                    k0,
                    a0,
                    b0,
                    ngr);

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

      if(any(group_log == g)){
        mujoin(g) = arma::join_cols(mu(g), mutemp(g));
        s2join(g) = arma::join_cols(s2(g), s2temp(g));
        int nkpt = ptilde(g).n_elem - 1;
        arma::vec index = arma::regspace(0, nkpt);
        probjoin(g) = arma::join_cols(ptilde(g).elem(arma::find(index < nkpt)),
                 ptilde(g)(nkpt) * freqtemp(g) / napprox);
      } else {
        ptilde(g).resize(2);
        ptilde(g).fill(1.0);

        mujoin(g) = mutemp(g);
        s2join(g) = s2temp(g);
        probjoin(g) = freqtemp(g) / napprox;
      }

    }

    for(arma::uword g = 1; g <= ngr; g++){

      mujoin_complete(g - 1)   = arma::join_cols(mujoin(g), mujoin(0));
      s2join_complete(g - 1)   = arma::join_cols(s2join(g), s2join(0));
      probjoin_complete(g - 1) = arma::join_cols(w(g - 1) * probjoin(g), (1 - w(g - 1)) * probjoin(0));

      max_val(g - 1) = mujoin(g).n_elem;
    }


    // update cluster allocations
    clust_update_DDP(data,
                      group,
                      group_log,
                      mujoin_complete,
                      s2join_complete,
                      probjoin_complete,
                      clust,
                      temp_proc_cum,
                      max_val,
                      iter,
                      ngr);

    mu = mujoin;
    s2 = s2join;

    // if the burn-in phase is complete
    if(iter >= nburn){

      // save clust and group_log in the export objects
      result_clust.row(iter - nburn) = arma::trans(arma::conv_to<arma::vec>::from(clust));
      result_group_log.row(iter - nburn)  = arma::trans(arma::conv_to<arma::vec>::from(group_log));

      // compute the density for each urn
      for(arma::uword g = 0; g < ngr; g++){
        if(out_dens){

          dens.col(g) = eval_density(grid, mujoin_complete(g),
                   s2join_complete(g), probjoin_complete(g));
        }
      }

      // update the slice of the current iteration in the density object
      if(light_dens){
        result_dens.slice(0) += dens;
      } else {
        result_dens.slice(iter - nburn) = dens;
      }
      res_wei.row(iter - nburn) = w.t();

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
    Rcpp::Rcout << "\n" << "Estimation done in " << double(end_s-start_s)/CLOCKS_PER_SEC << " seconds\n";
  }

  // export a list with the results
  Rcpp::List resu;
  if(light_dens){
    resu["dens"]   = result_dens / (niter - nburn);
  } else {
    resu["dens"]   = result_dens;
  }
  resu["clust"]  = result_clust;
  resu["group_log"]   = result_group_log;
  resu["time"]   = double(end_s-start_s)/CLOCKS_PER_SEC;
  resu["wvals"]  = res_wei;

  return resu;
}
