#ifndef DDP_TWO_H
#define DDP_TWO_H

void accelerate_DDP2(arma::vec data,
                     arma::vec group,
                     arma::vec group_log,
                     arma::vec col_log,
                     arma::vec row_log,
                     arma::field<arma::vec> &mu,
                     arma::field<arma::vec> &s2,
                     arma::field<arma::vec> &muc,
                     arma::field<arma::vec> &s2c,
                     arma::field<arma::vec> &mur,
                     arma::field<arma::vec> &s2r,
                     arma::vec clust,
                     double m0,
                     double k0,
                     double a0,
                     double b0,
                     int ngr,
                     int ngrc,
                     int ngrr);

void para_clean_DDP2(arma::field<arma::vec> &mu,
                     arma::field<arma::vec> &s2,
                     arma::field<arma::vec> &muc,
                     arma::field<arma::vec> &s2c,
                     arma::field<arma::vec> &mur,
                     arma::field<arma::vec> &s2r,
                     arma::vec &clust,
                     arma::vec group,
                     arma::vec zeta,
                     arma::vec col_log,
                     arma::vec row_log,
                     int ngr,
                     int ngrc,
                     int ngrr);

void simu_trunc_DDP2(arma::field<arma::vec> &mutemp,
                     arma::field<arma::vec> &s2temp,
                     arma::field<arma::vec> &freqtemp,
                     arma::field<arma::vec> &mutempc,
                     arma::field<arma::vec> &s2tempc,
                     arma::field<arma::vec> &freqtempc,
                     arma::field<arma::vec> &mutempr,
                     arma::field<arma::vec> &s2tempr,
                     arma::field<arma::vec> &freqtempr,
                     double mass,
                     double wei_group,
                     double wei_col,
                     double m0,
                     double k0,
                     double a0,
                     double b0,
                     int napprox,
                     int ngr,
                     int ngrc,
                     int ngrr);

void clust_update_DDP2(arma::vec data,
                       arma::vec group,
                       arma::vec col_group,
                       arma::vec row_group,
                       arma::vec &group_log,
                       arma::vec &col_log,
                       arma::vec &row_log,
                       arma::field<arma::vec> mujoin_complete,
                       arma::field<arma::vec> s2join_complete,
                       arma::field<arma::vec> probjoin_complete,
                       arma::vec &clust,
                       arma::vec max_val_group,
                       arma::vec max_val_col,
                       int ngr,
                       int ngrc,
                       int ngrr);

void update_w_DDP2(arma::vec &w,
                   double mass,
                   double wei_group,
                   int n_approx_unif,
                   arma::vec group_log,
                   arma::vec clust,
                   arma::vec group,
                   int ngr);

void update_v_DDP2(arma::vec &v,
                   double mass,
                   double wei_group,
                   double wei_col,
                   int n_approx_unif,
                   arma::vec col_log,
                   arma::vec clust,
                   arma::vec group,
                   arma::vec col_group,
                   int ngr,
                   int ngrc);

#endif
