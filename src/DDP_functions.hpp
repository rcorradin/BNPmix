#ifndef DDP_H
#define DDP_H

void accelerate_DDP(arma::vec data,
                    arma::vec group,
                    arma::vec zeta,
                    arma::field<arma::vec> &mu,
                    arma::field<arma::vec> &s2,
                    arma::vec clust,
                    double m0,
                    double k0,
                    double a0,
                    double b0,
                    int ngr);

void para_clean_DDP(arma::field<arma::vec> &mu,
                    arma::field<arma::vec> &s2,
                    arma::vec &clust,
                    arma::vec group,
                    arma::vec zeta,
                    int ngr);

void simu_trunc_DDP(arma::field<arma::vec> &mutemp,
                    arma::field<arma::vec> &s2temp,
                    arma::field<arma::vec> &freqtemp,
                    double mass,
                    double wei,
                    double m0,
                    double k0,
                    double a0,
                    double b0,
                    int napprox,
                    int ngr);

void clust_update_DDP(arma::vec data,
                      arma::vec group,
                      arma::vec zeta,
                      arma::field<arma::vec> mujoin,
                      arma::field<arma::vec> s2join,
                      arma::field<arma::vec> probjoin,
                      arma::vec &clust,
                      arma::vec max_val,
                      int iter,
                      int ngr,
                      arma::vec &new_val);

void zeta_update_DDP(arma::vec &zeta,
                     arma::vec group,
                     arma::vec &clust,
                     arma::field<arma::vec> &mu,
                     arma::field<arma::vec> &s2,
                     arma::field<arma::vec> &ptilde,
                     arma::vec w,
                     double mass,
                     double wei,
                     int ngr);  

void update_w_DDP(arma::vec &w,
                  double mass,
                  double wei,
                  int n_approx_unif,
                  arma::vec zeta,
                  arma::vec clust,
                  arma::vec group,
                  int ngr);

#endif