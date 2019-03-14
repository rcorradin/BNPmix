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
#include "IcsFunctions.h"
// [[Rcpp::depends("RcppArmadillo")]]


//' @export cICS
//' @name cICS
//' @title C++ function to estimate Pitman-Yor univariate mixtures via importance conditional sampler
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
//' @param napprox number of approximating values
//' @param nupd number of iterations to show current updating
//' @param out_param if TRUE, return also the location and scale paramteres lists
//' @param out_dens if TRUE, return also the estimated density (default TRUE)
//' @param sigma_PY second parameter of PY
//' @param print_message print the status

//[[Rcpp::export]]
Rcpp::List cICS(arma::vec data,
                arma::vec grid,
                int niter,
                int nburn,
                double m0,
                double k0,
                double a0,
                double b0,
                double mass,
                int napprox,
                int nupd = 0,
                bool out_param = 0,
                bool out_dens = 1,
                double sigma_PY = 0,
                bool print_message = 1){

  if(nupd == 0){
    nupd = (int) (niter / 10);
  }
  int n = data.n_elem;

  // initialize results objects
  arma::mat result_clust(niter - nburn, n);
  std::list<arma::vec> result_mu;
  std::list<arma::vec> result_s2;
  std::list<arma::vec> result_probs;
  arma::mat result_dens(niter - nburn, grid.n_elem);
  result_dens.fill(0);

  // initialize required object inside the loop
  arma::vec clust(n);
  arma::vec mu(1);
  arma::vec s2(1);
  arma::vec mutemp(1);
  arma::vec s2temp(1);
  arma::vec freqtemp(1);
  arma::vec ptilde;
  arma::vec dens(grid);
  arma::vec new_val(niter);
  arma::vec n_clust(niter - nburn);

  // fill the initialized quantity
  clust.fill(0);
  mu.fill(m0);
  s2.fill(b0 / (a0 - 1));
  mutemp.fill(0);
  s2temp.fill(0);
  freqtemp.fill(1);
  new_val.fill(0);
  n_clust.fill(0);

  // time quantities
  int start_s = clock();
  int current_s;

  // strarting loop
  for(arma::uword iter = 0; iter < niter; iter++){

    // acceleration step
    accelerate_ICS(data,
                   mu,
                   s2,
                   clust,
                   m0,
                   k0,
                   a0,
                   b0);

    // sample the probability vector
    // from Dirichlet distribution
    arma::vec temp_freq = freq_vec(clust);
    ptilde = rdirich_mass(temp_freq - sigma_PY, mass + sigma_PY * temp_freq.n_elem);

    // simulate the required values
    simu_trunc_PY(mutemp,
                  s2temp,
                  freqtemp,
                  mass + sigma_PY * temp_freq.n_elem,
                  m0,
                  k0,
                  a0,
                  b0,
                  napprox,
                  sigma_PY);

    // joint the existent parameters with
    // the simulated new ones
    arma::vec mujoin = arma::join_cols(mu, mutemp);
    arma::vec s2join = arma::join_cols(s2, s2temp);
    int nkpt = ptilde.n_elem;
    arma::vec index = arma::linspace(0, nkpt - 1, nkpt);
    arma::vec probjoin = arma::join_cols(ptilde.elem(arma::find(index < nkpt - 1)),
                                         ptilde[nkpt - 1] * freqtemp / napprox);

    // update cluster allocations
    clust_update_ICS(data,
                     mujoin,
                     s2join,
                     probjoin,
                     clust,
                     mu.n_elem,
                     iter,
                     new_val);

    mu = mujoin;
    s2 = s2join;

    // if the burn-in phase is complete
    if(iter >= nburn){
      result_mu.push_back(mu);
      result_s2.push_back(s2);
      result_probs.push_back(probjoin);
      n_clust(iter - nburn) = probjoin.n_elem;
      if(out_dens){
        dens = eval_density(grid,
                            mujoin,
                            s2join,
                            probjoin);
        result_dens.row(iter - nburn) = arma::trans(dens);
      }
    }

    // clean parameter objects
    para_clean_ICS(mu,
                   s2,
                   clust);

    // if the burn-in phase is complete
    if(iter >= nburn){
      result_clust.row(iter - nburn) = arma::trans(arma::conv_to<arma::vec>::from(clust));
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

//' @export cICS_mv
//' @name cICS_mv
//' @title C++ function to estimate Pitman-Yor multivariate mixtures via importance conditional sampler
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
//' @param napprox number of approximating values
//' @param nupd number of iterations to show current updating
//' @param out_param if TRUE, return also the location and scale paramteres lists
//' @param out_dens if TRUE, return also the estimated density (default TRUE)
//' @param sigma_PY second parameter of PY
//' @param print_message print the status
//' @param light_dens if TRUE return only the posterior mean of the density

//[[Rcpp::export]]
Rcpp::List cICS_mv(arma::mat data,
                   arma::mat grid,
                   int niter,
                   int nburn,
                   arma::vec m0,
                   double k0,
                   arma::mat S0,
                   double n0,
                   double mass,
                   int napprox,
                   int nupd = 0,
                   bool out_param = 0,
                   bool out_dens = 1,
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
  std::list<arma::mat> result_mu;
  std::list<arma::cube> result_s2;
  std::list<arma::vec> result_probs;
  arma::mat result_dens(grid.n_rows,1);
  if(!light_dens){
    result_dens.resize(niter - nburn, grid.n_rows);
  }
  result_dens.fill(0);
  arma::vec n_clust(niter - nburn);

  // initialize required object inside the loop
  arma::vec clust(n);
  arma::mat mu(1,d);
  arma::cube s2(d,d,1);
  arma::mat mutemp(1,d);
  arma::cube s2temp(d,d,1);
  arma::vec freqtemp(1);
  arma::vec ptilde;
  arma::vec dens(grid.n_rows);
  arma::vec new_val(niter);

  // fill the initialized quantity
  clust.fill(0);
  mu.row(0) = arma::trans(m0);
  s2.slice(0) = S0 / (n0 - d - 1);
  mutemp.fill(0);
  s2temp.fill(0);
  freqtemp.fill(1);
  new_val.fill(0);
  n_clust.fill(0);

  // time quantities
  int start_s = clock();
  int current_s;

  // strarting loop
  for(arma::uword iter = 0; iter < niter; iter++){

    // acceleration step
    accelerate_ICS_mv(data,
                      mu,
                      s2,
                      clust,
                      m0,
                      k0,
                      S0,
                      n0);

    // sample the probability vector
    // from Dirichlet distribution
    arma::vec temp_freq = freq_vec(clust);
    ptilde = rdirich_mass(temp_freq - sigma_PY, mass + sigma_PY * temp_freq.n_elem);

    // simulate the required values
    simu_trunc_PY_mv(mutemp,
                     s2temp,
                     freqtemp,
                     mass + sigma_PY * temp_freq.n_elem,
                     m0,
                     k0,
                     S0,
                     n0,
                     napprox,
                     sigma_PY);

    // joint the existent parameters with
    // the simulated new ones
    arma::mat mujoin = arma::join_cols(mu, mutemp);
    arma::cube s2join = arma::join_slices(s2, s2temp);
    int nkpt = ptilde.n_elem;
    arma::vec index = arma::linspace(0, nkpt - 1, nkpt);
    arma::vec probjoin = arma::join_cols(ptilde.elem(arma::find(index < nkpt - 1)),
                                         ptilde[nkpt - 1] * freqtemp / napprox);

    // update cluster allocations
    clust_update_ICS_mv(data,
                        mujoin,
                        s2join,
                        probjoin,
                        clust,
                        mu.n_elem,
                        iter,
                        new_val);

    mu = mujoin;
    s2 = s2join;

    // if the burn-in phase is complete
    if(iter >= nburn){
      result_mu.push_back(mu);
      result_s2.push_back(s2);
      result_probs.push_back(probjoin);
      n_clust(iter - nburn) = probjoin.n_elem;
      if(out_dens){
        dens = eval_density_mv(grid,
                               mujoin,
                               s2join,
                               probjoin);
        if(light_dens){
          result_dens += dens;
        } else {
          result_dens.row(iter - nburn) = arma::trans(dens);
        }
      }
    }

    // clean parameter objects
    para_clean_ICS_mv(mu,
                      s2,
                      clust);

    // if the burn-in phase is complete
    if(iter >= nburn){
      result_clust.row(iter - nburn) = arma::trans(arma::conv_to<arma::vec>::from(clust));
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
