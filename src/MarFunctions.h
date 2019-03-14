#ifndef MAR_H
#define MAR_H

void accelerate_MAR(arma::vec data,
                    arma::vec &mu,
                    arma::vec &s2,
                    arma::vec clust,
                    double m0,
                    double k0,
                    double a0,
                    double b0);

void accelerate_MAR_mv(arma::mat data,
                       arma::mat &mu,
                       arma::cube &s2,
                       arma::vec clust,
                       arma::vec m0,
                       double k0,
                       arma::mat S0,
                       double n0);

void para_clean_MAR(arma::vec &mu,
                    arma::vec &s2,
                    arma::vec &clust);

void para_clean_MAR_mv(arma::mat &mu,
                       arma::cube &s2,
                       arma::vec &clust);

void clust_update_MAR_PY(arma::vec data,
                         arma::vec &mu,
                         arma::vec &s2,
                         arma::vec &clust,
                         double mass,
                         double m0,
                         double k0,
                         double a0,
                         double b0,
                         int iter,
                         arma::vec &new_val,
                         double sigma_PY);

void clust_update_MAR_PY_mv(arma::mat data,
                            arma::mat &mu,
                            arma::cube &s2,
                            arma::vec &clust,
                            double mass,
                            arma::vec m0,
                            double k0,
                            arma::mat S0,
                            double n0,
                            int iter,
                            arma::vec &new_val,
                            double sigma_PY);

#endif
