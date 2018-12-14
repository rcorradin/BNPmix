#ifndef cSLI_H
#define cSLI_H

void accelerate_SLI(arma::vec data,
                    arma::vec &mu,
                    arma::vec &s2,
                    arma::vec &v,
                    arma::vec &w,
                    arma::vec clust,
                    double m0,
                    double k0,
                    double a0,
                    double b0,
                    double mass);

void accelerate_SLI_PY(arma::vec data,
                       arma::vec &mu,
                       arma::vec &s2,
                       arma::vec &v,
                       arma::vec &w,
                       arma::vec clust,
                       double m0,
                       double k0,
                       double a0,
                       double b0,
                       double mass,
                       double sigma_PY);

void accelerate_SLI_mv(arma::mat data,
                       arma::mat &mu,
                       arma::cube &s2,
                       arma::vec &v,
                       arma::vec &w,
                       arma::vec clust,
                       arma::vec m0,
                       double k0,
                       arma::mat S0,
                       double n0,
                       double mass);

void accelerate_SLI_PY_mv(arma::mat data,
                          arma::mat &mu,
                          arma::cube &s2,
                          arma::vec &v,
                          arma::vec &w,
                          arma::vec clust,
                          arma::vec m0,
                          double k0,
                          arma::mat S0,
                          double n0,
                          double mass,
                          double sigma_PY);

void para_clean_SLI(arma::vec &mu,
                    arma::vec &s2,
                    arma::vec &clust,
                    arma::vec &v,
                    arma::vec &w);

void para_clean_SLI_mv(arma::mat &mu,
                       arma::cube &s2,
                       arma::vec &clust,
                       arma::vec &v,
                       arma::vec &w);
                      
void grow_param_SLI(arma::vec &mu,
                arma::vec &s2,
                arma::vec &v,
                arma::vec &w,
                arma::vec u,
                double m0,
                double k0,
                double a0,
                double b0,
                double mass,
                int n);

void grow_param_SLI_PY(arma::vec &mu,
                       arma::vec &s2,
                       arma::vec &v,
                       arma::vec &w,
                       arma::vec u,
                       double m0,
                       double k0,
                       double a0,
                       double b0,
                       double mass,
                       int n,
                       double sigma_PY);

void grow_param_SLI_mv(arma::mat &mu,
                       arma::cube &s2,
                       arma::vec &v,
                       arma::vec &w,
                       arma::vec u,
                       arma::vec m0,
                       double k0,
                       arma::mat S0,
                       double n0,
                       double mass,
                       int n);

void grow_param_SLI_PY_mv(arma::mat &mu,
                          arma::cube &s2,
                          arma::vec &v,
                          arma::vec &w,
                          arma::vec u,
                          arma::vec m0,
                          double k0,
                          arma::mat S0,
                          double n0,
                          double mass,
                          int n,
                          double sigma_PY);
                
void update_u_SLI(arma::vec clust,
              arma::vec w,
              arma::vec &u);
              
void update_cluster_SLI(arma::vec data,
                        arma::vec mu,
                        arma::vec s2,
                        arma::vec &clust,
                        arma::vec w,
                        arma::vec u,
                        int max_val,
                        int iter,
                        arma::vec &new_val);

void update_cluster_SLI_mv(arma::mat data,
                           arma::mat mu,
                           arma::cube s2,
                           arma::vec &clust,
                           arma::vec w,
                           arma::vec u,
                           int max_val,
                           int iter,
                           arma::vec &new_val);
                        
#endif