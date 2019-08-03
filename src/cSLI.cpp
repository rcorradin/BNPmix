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
#include "Distributions.h"
#include "CommonUtilities.h"
#include "SliFunctions.h"
// [[Rcpp::depends("RcppArmadillo")]]

//' @export
//' @name cSLI_L
//' @title C++ function to estimate Pitman-Yor univariate mixtures via slice sampler - LOCATION
//' @keywords internal
//'
//' @param data a vector of observations
//' @param grid vector to evaluate the density
//' @param niter number of iterations
//' @param nburn number of burn-in iterations
//' @param m0 expectation of location component
//' @param s20 variance of location component
//' @param a0 parameter of scale component
//' @param b0 parameter of scale component
//' @param m1 hyperparameter, mean of distribution of m0
//' @param k1 hyperparameter, scale factor of distribution of m0
//' @param a1 hyperparameter, shape of distribution of s20
//' @param b1 hyperparameter, rate of distribution of s20
//' @param mass mass parameter
//' @param nupd number of iterations to show current updating
//' @param out_param if TRUE, return also the location and scale paramteres lists
//' @param out_dens if TRUE, return also the estimated density (default TRUE)
//' @param sigma_PY second parameter of PY
//' @param print_message print the status
//' @param hyper, if TRUE use hyperpriors, default TRUE
//

//[[Rcpp::export]]
Rcpp::List cSLI_L(arma::vec data,
                  arma::vec grid,
                  int niter,
                  int nburn,
                  double m0,
                  double s20,
                  double a0,
                  double b0,
                  double m1,
                  double k1,
                  double a1,
                  double b1,
                  double mass,
                  int nupd = 0,
                  bool out_param = 0,
                  bool out_dens = 1,
                  double sigma_PY = 0,
                  bool print_message = 1,
                  bool hyper = true){

  if(nupd == 0){
    nupd = int(niter / 10);
  }

  int n = data.n_elem;

  arma::mat result_clust(niter - nburn, n);
  std::list<arma::vec> result_mu;
  arma::vec result_s2(niter - nburn);
  std::list<arma::vec> result_probs;
  arma::mat result_dens(niter - nburn, grid.n_elem);
  arma::vec new_val(niter);
  arma::vec n_clust(niter - nburn);
  result_dens.fill(0);

  arma::vec clust(n);
  clust.fill(0);

  arma::vec mu(1);
  double s2;
  arma::vec v(1);
  arma::vec w(1);
  arma::vec u(n, arma::fill::randu);

  new_val.fill(0);
  n_clust.fill(0);
  mu.fill(m0);
  s2 = (b0 / (a0 - 1));

  // time quantities
  int start_s = clock();
  int current_s;

  // strarting loop
  for(arma::uword iter = 0; iter < niter; iter++){

    // update the parameters
    accelerate_SLI_PY_L(data,
                        mu,
                        s2,
                        v,
                        w,
                        clust,
                        m0,
                        s20,
                        a0,
                        b0,
                        mass,
                        sigma_PY);

    // hyper-accelerate
    if(hyper){
      hyper_accelerate_SLI_L(mu,
                             clust,
                             m0,
                             s20,
                             m1,
                             k1,
                             a1,
                             b1);
    }

    int old_length = mu.n_elem;

    // update the slice weights
    update_u_SLI(clust,
                 w,
                 u);

    // extend the stick breaking representation
    grow_param_SLI_PY_L(mu,
                        v,
                        w,
                        u,
                        m0,
                        s20,
                        mass,
                        n,
                        sigma_PY);

    // update the allocation
    update_cluster_SLI_L(data,
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
      result_s2(iter - nburn) = s2;
      result_probs.push_back(w);
      n_clust(iter - nburn) = w.n_elem;
      if(out_dens){
        result_dens.row(iter - nburn) = arma::trans(eval_density_L(grid, mu, s2, w));
      }
    }

    mu.resize(max(clust) + 1);
    w.resize(max(clust) + 1);
    v.resize(max(clust) + 1);

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
    resu["time"]   = double(end_s-start_s)/CLOCKS_PER_SEC;
  } else {
    resu["dens"]   = result_dens;
    resu["clust"]  = result_clust;
    resu["time"]   = double(end_s-start_s)/CLOCKS_PER_SEC;
  }
  return resu;
}

//' @export
//' @name cSLI
//' @title C++ function to estimate Pitman-Yor univariate mixtures via slice sampler - LOCATION SCALE
//' @keywords internal
//'
//' @param data a vector of observations
//' @param grid vector to evaluate the density
//' @param niter number of iterations
//' @param nburn number of burn-in iterations
//' @param m0 expectation of location component
//' @param k0 tuning parameter of variance of location component
//' @param a0 parameter of scale component
//' @param b0 parameter of scale component
//' @param m1 mean of hyperdistribution of m0
//' @param s21 variance of hyperdistribution of m0
//' @param tau1 shape parameter of hyperdistribution of k0
//' @param tau2 rate parameter of hyperdistribution of k0
//' @param a1 shape parameter of hyperdistribution of b0
//' @param b1 rate parameter of hyperdistribution of b0
//' @param mass parameter
//' @param nupd number of iterations to show current updating
//' @param out_param if TRUE, return also the location and scale paramteres lists
//' @param out_dens if TRUE, return also the estimated density (default TRUE)
//' @param sigma_PY second parameter of PY
//' @param print_message print the status
//' @param hyper, if TRUE use hyperpriors, default TRUE
//

//[[Rcpp::export]]
Rcpp::List cSLI(arma::vec data,
                arma::vec grid,
                int niter,
                int nburn,
                double m0,
                double k0,
                double a0,
                double b0,
                double m1,
                double s21,
                double tau1,
                double tau2,
                double a1,
                double b1,
                double mass,
                int nupd = 0,
                bool out_param = 0,
                bool out_dens = 1,
                double sigma_PY = 0,
                bool print_message = 1,
                bool hyper = 1){
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

    if(hyper){
      hyper_accelerate_SLI_L(mu,
                             s2,
                             clust,
                             m0,
                             k0,
                             a0,
                             b0,
                             m1,
                             s21,
                             tau1,
                             tau2,
                             a1,
                             b1);
    }

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
    // para_clean_SLI(mu,
    //                s2,
    //                clust,
    //                v,
    //                w);


    mu.resize(max(clust) + 1);
    s2.resize(max(clust) + 1);
    w.resize(max(clust) + 1);
    v.resize(max(clust) + 1);

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
    resu["time"]   = double(end_s-start_s)/CLOCKS_PER_SEC;
  } else {
    resu["dens"]   = result_dens;
    resu["clust"]  = result_clust;
    resu["time"]   = double(end_s-start_s)/CLOCKS_PER_SEC;
  }
  return resu;
}

/*----------------------------------------------------------------------
 *
 * MULTIVARIATE
 * LOCATION KERNEL
 * SLI main
 *
 *----------------------------------------------------------------------
 */

//' @export
//' @name cSLI_mv_L
//' @title C++ function to estimate Pitman-Yor multivariate mixtures via slice sampler - LOCATION
//' @keywords internal
//'
//' @param data a matrix of observations
//' @param grid matrix of points to evaluate the density
//' @param niter number of iterations
//' @param nburn number of burn-in iterations
//' @param m0 expectation of location component
//' @param S20 variance of location component
//' @param S0 parameter of scale component
//' @param n0 parameter of scale component
//' @param m1 mean of hyperdistribtion of m0
//' @param k1 scale factor of hyperdistribtion of m0
//' @param theta1 df of hyperdistribtion of S20
//' @param Theta1 matrix of hyperdistribution of S20
//' @param mass mass parameter
//' @param nupd number of iterations to show current updating
//' @param out_param if TRUE, return also the location and scale paramteres lists
//' @param out_dens if TRUE, return also the estimated density (default TRUE)
//' @param sigma_PY second parameter of PY
//' @param print_message print the status
//' @param light_dens if TRUE return only the posterior mean of the density
//' @param hyper, if TRUE use hyperpriors, default TRUE
//

//[[Rcpp::export]]
Rcpp::List cSLI_mv_L(arma::mat data,
                     arma::mat grid,
                     int niter,
                     int nburn,
                     arma::vec m0,
                     arma::mat S20,
                     arma::mat S0,
                     double n0,
                     arma::vec m1,
                     double k1,
                     double theta1,
                     arma::mat Theta1,
                     double mass,
                     int nupd = 0,
                     bool out_param = 0,
                     bool out_dens = 1,
                     double sigma_PY = 0,
                     bool print_message = 1,
                     bool light_dens = 1,
                     bool hyper = 1){

  if(nupd == 0){
    nupd = int(niter / 10);
  }

  int n = data.n_rows;
  int d = data.n_cols;

  // initialize results objects
  arma::mat result_clust(niter - nburn, n);
  std::list<arma::mat> result_mu;
  arma::cube result_s2(d, d, niter - nburn);
  std::list<arma::vec> result_probs;
  arma::mat result_dens(grid.n_rows, 1);
  if(!light_dens){
    result_dens.resize(niter - nburn, grid.n_rows);
  }
  result_dens.fill(0);

  // initialize required object inside the loop
  arma::vec clust(n);
  arma::mat mu(1,d);
  arma::mat s2(d,d);
  arma::vec v(1);
  arma::vec w(1);
  arma::vec u(n, arma::fill::randu);
  arma::vec dens(grid.n_rows);
  arma::vec new_val(niter);
  arma::vec n_clust(niter - nburn);

  // fill the initialized quantity
  clust.fill(0);
  mu.row(0) = arma::trans(m0);
  s2 = S0 / (n0 - d - 1);
  new_val.fill(0);
  n_clust.fill(0);

  // time quantities
  int start_s = clock();
  int current_s;

  // strarting loop
  for(arma::uword iter = 0; iter < niter; iter++){

    // update the parameters
    accelerate_SLI_PY_mv_L(data,
                           mu,
                           s2,
                           v,
                           w,
                           clust,
                           m0,
                           S20,
                           S0,
                           n0,
                           mass,
                           sigma_PY);

    if(hyper){
      hyper_accelerate_MAR_mv_L(mu,
                                m0,
                                clust,
                                S20,
                                m1,
                                k1,
                                theta1,
                                Theta1);
    }

    int old_length = mu.n_elem;

    // update the stick breaking weights
    update_u_SLI(clust,
                 w,
                 u);

    // extend the stick breaking representation
    grow_param_SLI_PY_mv_L(mu,
                           v,
                           w,
                           u,
                           m0,
                           S20,
                           mass,
                           n,
                           sigma_PY);

    // update the allocation
    update_cluster_SLI_mv_L(data,
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
      result_s2.slice(iter - nburn) = (s2);
      result_probs.push_back(w);
      n_clust(iter - nburn) = w.n_elem;
      if(out_dens){
        dens = eval_density_mv_L(grid,
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

    mu.resize(max(clust) + 1, d);
    w.resize(max(clust) + 1);
    v.resize(max(clust) + 1);

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
    resu["time"]   = double(end_s-start_s)/CLOCKS_PER_SEC;
  } else {
    if(light_dens){
      resu["dens"]   = result_dens / (niter - nburn);
    } else {
      resu["dens"]   = result_dens;
    }
    resu["clust"]  = result_clust;
    resu["time"]   = double(end_s-start_s)/CLOCKS_PER_SEC;
  }
  return resu;
}

/*----------------------------------------------------------------------
 *
 * MULTIVARIATE
 * LOCATION-SCALE KERNEL
 * slice sampler functions
 *
 *----------------------------------------------------------------------
 */

//' @export
//' @name cSLI_mv
//' @title C++ function to estimate Pitman-Yor multivariate mixtures via slice sampler - LOCATION SCALE
//' @keywords internal
//'
//' @param data a matrix of observations
//' @param grid matrix of points to evaluate the density
//' @param niter number of iterations
//' @param nburn number of burn-in iterations
//' @param m0 expectation of location component
//' @param k0 tuning parameter of variance of location component
//' @param S0 parameter of scale component
//' @param n0 parameter of scale component
//' @param m1 mean of hyperprior distribution of m0
//' @param S1 covariance of hyperprior distribution of m0
//' @param tau1 shape parameter of hyperprior distribution of k0
//' @param tau2 rate parameter of hyperprior distribution of k0
//' @param theta1 df of hyperprior distribution of S0
//' @param Theta1 matrix of hyperprior distribution of S0
//' @param mass mass parameter
//' @param nupd number of iterations to show current updating
//' @param out_param if TRUE, return also the location and scale paramteres lists
//' @param out_dens if TRUE, return also the estimated density (default TRUE)
//' @param sigma_PY second parameter of PY
//' @param print_message print the status
//' @param light_dens if TRUE return only the posterior mean of the density
//' @param hyper, if TRUE use hyperpriors, default TRUE
//

//[[Rcpp::export]]
Rcpp::List cSLI_mv(arma::mat data,
                   arma::mat grid,
                   int niter,
                   int nburn,
                   arma::vec m0,
                   double k0,
                   arma::mat S0,
                   double n0,
                   arma::vec m1,
                   arma::mat S1,
                   double tau1,
                   double tau2,
                   double theta1,
                   arma::mat Theta1,
                   double mass,
                   int nupd = 0,
                   bool out_param = 0,
                   bool out_dens = 1,
                   double sigma_PY = 0,
                   bool print_message = 1,
                   bool light_dens = 1,
                   bool hyper = 1){

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

    if(hyper){
      hyper_accelerate_SLI_mv_LS(mu,
                                 s2,
                                 clust,
                                 m0,
                                 k0,
                                 S0,
                                 n0,
                                 m1,
                                 S1,
                                 tau1,
                                 tau2,
                                 theta1,
                                 Theta1);
    }

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

    mu.resize(max(clust) + 1, d);
    s2.resize(d, d, max(clust) + 1);
    w.resize(max(clust) + 1);
    v.resize(max(clust) + 1);

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
    resu["time"]   = double(end_s-start_s)/CLOCKS_PER_SEC;
  } else {
    if(light_dens){
      resu["dens"]   = result_dens / (niter - nburn);
    } else {
      resu["dens"]   = result_dens;
    }
    resu["clust"]  = result_clust;
    resu["time"]   = double(end_s-start_s)/CLOCKS_PER_SEC;
  }
  return resu;
}

/*----------------------------------------------------------------------
 *
 * MULTIVARIATE
 * PRODUCT KERNEL
 * slice sampler functions
 *
 *----------------------------------------------------------------------
 */

//' @export
//' @name cSLI_mv_P
//' @title C++ function to estimate Pitman-Yor multivariate mixtures via slice sampler - LOCATION SCALE
//' @keywords internal
//'
//' @param data a matrix of observations
//' @param grid matrix of points to evaluate the density
//' @param niter number of iterations
//' @param nburn number of burn-in iterations
//' @param m0 expectation of location component
//' @param k0 tuning parameter of variance of location component
//' @param S0 parameter of scale component
//' @param n0 parameter of scale component
//' @param m1 means of hyperdistribution of m0
//' @param s21 variances of hyperdistribution of m0
//' @param tau1 shape parameters of hyperdistribution of k0
//' @param tau2 rate parameters of hyperdistribution of k0
//' @param a1 shape parameters of hyperdistribution of b0
//' @param b1 rate parameters of hyperdistribution of b0
//' @param mass mass parameter
//' @param nupd number of iterations to show current updating
//' @param out_param if TRUE, return also the location and scale paramteres lists
//' @param out_dens if TRUE, return also the estimated density (default TRUE)
//' @param sigma_PY second parameter of PY
//' @param print_message print the status
//' @param light_dens if TRUE return only the posterior mean of the density
//' @param hyper, if TRUE use hyperpriors, default TRUE
//

//[[Rcpp::export]]
Rcpp::List cSLI_mv_P(arma::mat data,
                     arma::mat grid,
                     int niter,
                     int nburn,
                     arma::vec m0,
                     arma::vec k0,
                     arma::vec a0,
                     arma::vec b0,
                     arma::vec m1,
                     arma::vec s21,
                     arma::vec tau1,
                     arma::vec tau2,
                     arma::vec a1,
                     arma::vec b1,
                     double mass,
                     int nupd = 0,
                     bool out_param = 0,
                     bool out_dens = 1,
                     double sigma_PY = 0,
                     bool print_message = 1,
                     bool light_dens = 1,
                     bool hyper = 1){

  if(nupd == 0){
    nupd = int(niter / 10);
  }

  int n = data.n_rows;
  int d = data.n_cols;

  // initialize results objects
  arma::mat result_clust(niter - nburn, n);
  std::list<arma::mat> result_mu;
  std::list<arma::mat> result_s2;
  std::list<arma::vec> result_probs;
  arma::mat result_dens(grid.n_rows, 1);
  if(!light_dens){
    result_dens.resize(niter - nburn, grid.n_rows);
  }
  result_dens.fill(0);

  // initialize required object inside the loop
  arma::vec clust(n);
  arma::mat mu(1,d);
  arma::mat s2(1,d);
  arma::vec v(1);
  arma::vec w(1);
  arma::vec u(n, arma::fill::randu);
  arma::vec dens(grid.n_rows);
  arma::vec new_val(niter);
  arma::vec n_clust(niter - nburn);

  // fill the initialized quantity
  clust.fill(0);
  mu.row(0) = arma::trans(m0);
  s2.row(0) = arma::trans(b0 / (a0 - 1));
  new_val.fill(0);
  n_clust.fill(0);

  // time quantities
  int start_s = clock();
  int current_s;

  // strarting loop
  for(arma::uword iter = 0; iter < niter; iter++){

    // update the parameters
    accelerate_SLI_PY_mv_P(data,
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

    if(hyper){
      hyper_accelerate_SLI_mv_P(mu,
                                s2,
                                clust,
                                m0,
                                k0,
                                a0,
                                b0,
                                m1,
                                s21,
                                tau1,
                                tau2,
                                a1,
                                b1);
    }

    int old_length = mu.n_elem;

    // update the stick breaking weights
    update_u_SLI(clust,
                 w,
                 u);

    // extend the stick breaking representation
    grow_param_SLI_PY_mv_P(mu,
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
    update_cluster_SLI_mv_P(data,
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
        dens = eval_density_mv_P(grid,
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

    mu.resize(max(clust) + 1, d);
    s2.resize(max(clust) + 1,d);
    w.resize(max(clust) + 1);
    v.resize(max(clust) + 1);

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
    resu["time"]   = double(end_s-start_s)/CLOCKS_PER_SEC;
  } else {
    if(light_dens){
      resu["dens"]   = result_dens / (niter - nburn);
    } else {
      resu["dens"]   = result_dens;
    }
    resu["clust"]  = result_clust;
    resu["time"]   = double(end_s-start_s)/CLOCKS_PER_SEC;
  }
  return resu;
}

/*----------------------------------------------------------------------
 *
 * MULTIVARIATE
 * MIXTURE OF KERNEL REGRESSION
 * slice sampler function
 *
 *----------------------------------------------------------------------
 */

//' @export
//' @name cSLI_mv_MKR
//' @title C++ function to estimate Pitman-Yor multivariate mixtures via slice sampler - PRODUCT KERNEL
//' @keywords internal
//'
//' @param data a matrix of observations
//' @param grid matrix of points to evaluate the density
//' @param niter number of iterations
//' @param nburn number of burn-in iterations
//' @param m0 expectation of location component
//' @param k0 vector, scale parameters for the location component
//' @param a0 vector, parameters of scale component
//' @param b0 vector, parameters of scale component
//' @param m1 means of hyperdistribution of m0
//' @param s21 variances of hyperdistribution of m0
//' @param tau1 shape parameters of hyperdistribution of k0
//' @param tau2 rate parameters of hyperdistribution of k0
//' @param a1 shape parameters of hyperdistribution of b0
//' @param b1 rate parameters of hyperdistribution of b0
//' @param strength strength parameter
//' @param napprox number of approximating values
//' @param nupd number of iterations to show current updating
//' @param out_param if TRUE, return also the location and scale paramteres lists
//' @param out_dens if TRUE, return also the estimated density (default TRUE)
//' @param discount second parameter of PY
//' @param print_message print the status
//' @param light_dens if TRUE return only the posterior mean of the density
//' @param hyper, if TRUE use hyperpriors, default TRUE
//

//[[Rcpp::export]]
Rcpp::List cSLI_mv_MKR(arma::vec y,
                       arma::mat covs,
                       arma::vec grid_response,
                       arma::mat grid_covs,
                       int niter,
                       int nburn,
                       arma::vec beta0,
                       arma::mat Sb0,
                       double a0,
                       double b0,
                       arma::vec beta1,
                       double k1,
                       double sb1,
                       arma::mat Sb1,
                       double tau1,
                       double tau2,
                       double strength,
                       int nupd = 0,
                       bool out_param = 0,
                       bool out_dens = 1,
                       double discount = 0,
                       bool print_message = 1,
                       bool light_dens = 1,
                       bool hyper = 1){
  if(nupd == 0){
    nupd = (int) (niter / 10);
  }
  int n = covs.n_rows;
  int d = covs.n_cols;

  // initialize results objects
  arma::mat result_clust(niter - nburn, n);
  arma::mat result_beta(niter - nburn, d);
  std::list<arma::mat> result_beta_all;
  std::list<arma::vec> result_sigma2;
  std::list<arma::vec> result_probs;

  // initialize the result densities objects (response and covs)
  arma::cube result_dens(grid_response.n_rows, grid_covs.n_rows, 1);
  if(!light_dens){
    result_dens.resize(grid_response.n_rows, grid_covs.n_rows, niter - nburn);
  }
  result_dens.fill(0);

  // initialize required object inside the loop
  arma::vec clust(n);
  arma::mat beta(1,d);
  arma::vec sigma2(1);
  arma::mat dens(grid_response.n_rows, grid_covs.n_rows);
  arma::vec v(1);
  arma::vec w(1);
  arma::vec u(n, arma::fill::randu);

  // fill the initialized quantity
  clust.fill(0);
  beta.row(0) = arma::trans(beta0);
  sigma2.row(0) = b0 / (a0 - 1);

  // time quantities
  int start_s = clock();
  int current_s;

  // strarting loop
  for(arma::uword iter = 0; iter < niter; iter++){

    // acceleration step
    accelerate_SLI_mv_MRK(y,
                          covs,
                          beta,
                          sigma2,
                          v,
                          w,
                          clust,
                          beta0,
                          Sb0,
                          a0,
                          b0,
                          strength,
                          discount);

    if(hyper){
      hyper_accelerate_SLI_mv_MRK(y,
                                  covs,
                                  clust,
                                  beta,
                                  sigma2,
                                  beta0,
                                  Sb0,
                                  a0,
                                  b0,
                                  beta1,
                                  k1,
                                  sb1,
                                  Sb1,
                                  tau1,
                                  tau2);

    }

    // update the stick breaking weights
    update_u_SLI(clust,
                 w,
                 u);

    // extend the stick breaking representation
    grow_param_SLI_PY_mv_MRK(beta,
                             sigma2,
                             v,
                             w,
                             u,
                             beta0,
                             Sb0,
                             a0,
                             b0,
                             strength,
                             n,
                             discount);

    // update the allocation
    update_cluster_SLI_mv_MRK(y,
                              covs,
                              beta,
                              sigma2,
                              clust,
                              w,
                              u);

    // if the burn-in phase is complete
    if(iter >= nburn){
      result_probs.push_back(w);
      result_beta_all.push_back(beta);
      result_sigma2.push_back(sigma2);
      result_clust.row(iter - nburn)  = arma::trans(arma::conv_to<arma::vec>::from(clust));

      if(out_dens){
        dens = eval_density_mv_MKR(grid_covs,
                                   grid_response,
                                   beta,
                                   sigma2,
                                   w);

        if(light_dens){
          result_dens.slice(0) += dens;
        } else {
          result_dens.slice(iter - nburn) = dens;
        }
      }
    }

    beta.resize(max(clust) + 1, d);
    sigma2.resize(max(clust) + 1);
    w.resize(max(clust) + 1);
    v.resize(max(clust) + 1);


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
    resu["beta"]   = result_beta_all;
    resu["sigma2"] = result_sigma2;
    resu["probs"]  = result_probs;
    resu["time"]   = double(end_s-start_s)/CLOCKS_PER_SEC;
  } else {
    if(light_dens){
      resu["dens"]   = result_dens / (niter - nburn);
    } else {
      resu["dens"]   = result_dens;
    }
    resu["clust"]  = result_clust;
    resu["time"]   = double(end_s-start_s)/CLOCKS_PER_SEC;
  }
  return resu;
}
