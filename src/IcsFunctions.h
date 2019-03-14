#ifndef cICS_H
#define cICS_H

void accelerate_ICS(arma::vec data,
                    arma::vec &mu,
                    arma::vec &s2,
                    arma::vec clust,
                    double m0,
                    double k0,
                    double a0,
                    double b0);

void accelerate_ICS_mv(arma::mat data,
                       arma::mat &mu,
                       arma::cube &s2,
                       arma::vec clust,
                       arma::vec m0,
                       double k0,
                       arma::mat S0,
                       double n0);

void para_clean_ICS(arma::vec &mu,
                    arma::vec &s2,
                    arma::vec &clust);

void para_clean_ICS_mv(arma::mat &mu,
                       arma::cube &s2,
                       arma::vec &clust);

void simu_trunc_PY(arma::vec &mutemp,
                    arma::vec &s2temp,
                    arma::vec &freqtemp,
                    double mass,
                    double m0,
                    double k0,
                    double a0,
                    double b0,
                    int napprox,
                    double sigma_PY);

void simu_trunc_PY_mv(arma::mat &mutemp,
                       arma::cube &s2temp,
                       arma::vec &freqtemp,
                       double mass,
                       arma::vec m0,
                       double k0,
                       arma::mat S0,
                       double n0,
                       int napprox,
                       double sigma_PY);

void clust_update_ICS(arma::vec data,
                      arma::vec mujoin,
                      arma::vec s2join,
                      arma::vec probjoin,
                      arma::vec &clust,
                      int max_val,
                      int iter,
                      arma::vec &new_val);

void clust_update_ICS_mv(arma::mat data,
                         arma::mat mujoin,
                         arma::cube s2join,
                         arma::vec probjoin,
                         arma::vec &clust,
                         int max_val,
                         int iter,
                         arma::vec &new_val);

#endif
