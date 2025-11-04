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
                      arma::vec &zeta,
                      arma::field<arma::vec> mujoin_complete,
                      arma::field<arma::vec> s2join_complete,
                      arma::field<arma::vec> probjoin_complete,
                      arma::vec &clust,
                      arma::mat &temp_proc_cum,
                      arma::vec max_val,
                      int iter,
                      int ngr);

void update_w_DDP(arma::vec &w,
                  double mass,
                  double wei,
                  double var_MH_step,
                  arma::vec zeta,
                  arma::vec clust,
                  arma::vec group,
                  arma::mat temp_proc_cum,
                  int ngr);

#endif
