#ifndef cSLI_H
#define cSLI_H

void update_u_SLI(arma::vec clust,
                  arma::vec w,
                  arma::vec &u);

/*----------------------------------------------------------------------
 *
 * UNIVARIATE
 * LOCATION KERNEL
 * ICS functions
 *
 *----------------------------------------------------------------------
 */

void accelerate_SLI_PY_L(arma::vec data,
                         arma::vec &mu,
                         double &s2,
                         arma::vec &v,
                         arma::vec &w,
                         arma::vec clust,
                         double m0,
                         double s20,
                         double a0,
                         double b0,
                         double mass,
                         double sigma_PY);

void hyper_accelerate_SLI_L(arma::vec mu,
                            arma::vec clust,
                            double &m0,
                            double &s20,
                            double m1,
                            double k1,
                            double a1,
                            double b1);

void grow_param_SLI_PY_L(arma::vec &mu,
                         arma::vec &v,
                         arma::vec &w,
                         arma::vec u,
                         double m0,
                         double s20,
                         double mass,
                         int n,
                         double sigma_PY);

void grow_param_indep_SLI_PY_L(arma::vec &mu,
                               arma::vec &v,
                               arma::vec &w,
                               arma::vec &xi,
                               arma::vec u,
                               double m0,
                               double s20,
                               double mass,
                               int n,
                               double sigma_PY,
                               double param_seq_one,
                               double param_seq_two);

void update_cluster_SLI_L(arma::vec data,
                          arma::vec mu,
                          double s2,
                          arma::vec &clust,
                          arma::vec w,
                          arma::vec u);


void update_cluster_indep_SLI_L(arma::vec data,
                                arma::vec mu,
                                double s2,
                                arma::vec &clust,
                                arma::vec w,
                                arma::vec xi,
                                arma::vec u);

/*----------------------------------------------------------------------
 *
 * UNIVARIATE
 * LOCATION-SCALE KERNEL
 * ICS functions
 *
 *----------------------------------------------------------------------
 */

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

void hyper_accelerate_SLI(arma::vec mu,
                          arma::vec s2,
                          arma::vec clust,
                          double &m0,
                          double &k0,
                          double a0,
                          double &b0,
                          double m1,
                          double s21,
                          double tau1,
                          double tau2,
                          double a1,
                          double b1);

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
                       double sigma_PY,
                       int &bound);

void grow_param_indep_SLI_PY(arma::vec &mu,
                             arma::vec &s2,
                             arma::vec &v,
                             arma::vec &w,
                             arma::vec &xi,
                             arma::vec u,
                             double m0,
                             double k0,
                             double a0,
                             double b0,
                             double mass,
                             int n,
                             double sigma_PY,
                             double param_seq_one,
                             double param_seq_two,
                             int &bound);

void update_cluster_SLI(arma::vec data,
                        arma::vec mu,
                        arma::vec s2,
                        arma::vec &clust,
                        arma::vec w,
                        arma::vec u);

void update_cluster_indep_SLI(arma::vec data,
                              arma::vec mu,
                              arma::vec s2,
                              arma::vec &clust,
                              arma::vec w,
                              arma::vec xi,
                              arma::vec u);

/*----------------------------------------------------------------------
 *
 * MULTIVARIATE
 * LOCATION KERNEL
 * ICS functions
 *
 *----------------------------------------------------------------------
 */

void accelerate_SLI_PY_mv_L(arma::mat data,
                            arma::mat &mu,
                            arma::mat &s2,
                            arma::vec &v,
                            arma::vec &w,
                            arma::vec clust,
                            arma::vec m0,
                            arma::mat S20,
                            arma::mat S0,
                            double n0,
                            double mass,
                            double sigma_PY);

void hyper_accelerate_SLI_mv_L(arma::mat mu,
                               arma::vec &m0,
                               arma::vec clust,
                               arma::mat &S20,
                               arma::vec m1,
                               double k1,
                               double theta1,
                               arma::mat Theta1);

void grow_param_SLI_PY_mv_L(arma::mat &mu,
                            arma::vec &v,
                            arma::vec &w,
                            arma::vec u,
                            arma::vec m0,
                            arma::mat S20,
                            double mass,
                            int n,
                            double sigma_PY);

void grow_param_indep_SLI_PY_mv_L(arma::mat &mu,
                                  arma::vec &v,
                                  arma::vec &w,
                                  arma::vec &xi,
                                  arma::vec u,
                                  arma::vec m0,
                                  arma::mat S20,
                                  double mass,
                                  int n,
                                  double sigma_PY,
                                  double param_seq_one,
                                  double param_seq_two);

void update_cluster_SLI_mv_L(arma::mat data,
                             arma::mat mu,
                             arma::mat s2,
                             arma::vec &clust,
                             arma::vec w,
                             arma::vec u);

void update_cluster_indep_SLI_mv_L(arma::mat data,
                                   arma::mat mu,
                                   arma::mat s2,
                                   arma::vec &clust,
                                   arma::vec w,
                                   arma::vec xi,
                                   arma::vec u);

/*----------------------------------------------------------------------
 *
 * MULTIVARIATE
 * LOCATION-SCALE KERNEL
 * ICS functions
 *
 *----------------------------------------------------------------------
 */

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

void hyper_accelerate_SLI_mv_LS(arma::mat mu,
                                arma::cube s2,
                                arma::vec clust,
                                arma::vec &m0,
                                double &k0,
                                arma::mat &S0,
                                double n0,
                                arma::vec m1,
                                arma::mat S1,
                                double tau1,
                                double tau2,
                                double theta1,
                                arma::mat Theta1);

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

void grow_param_indep_SLI_PY_mv(arma::mat &mu,
                                arma::cube &s2,
                                arma::vec &v,
                                arma::vec &w,
                                arma::vec &xi,
                                arma::vec u,
                                arma::vec m0,
                                double k0,
                                arma::mat S0,
                                double n0,
                                double mass,
                                int n,
                                double sigma_PY,
                                double param_seq_one,
                                double param_seq_two);

void update_cluster_SLI_mv(arma::mat data,
                           arma::mat mu,
                           arma::cube s2,
                           arma::vec &clust,
                           arma::vec w,
                           arma::vec u);

void update_cluster_indep_SLI_mv(arma::mat data,
                                 arma::mat mu,
                                 arma::cube s2,
                                 arma::vec &clust,
                                 arma::vec w,
                                 arma::vec xi,
                                 arma::vec u);

/*----------------------------------------------------------------------
 *
 * MULTIVARIATE
 * PRODUCT KERNEL
 * ICS functions
 *
 *----------------------------------------------------------------------
 */

void accelerate_SLI_PY_mv_P(arma::mat data,
                            arma::mat &mu,
                            arma::mat &s2,
                            arma::vec &v,
                            arma::vec &w,
                            arma::vec clust,
                            arma::vec m0,
                            arma::vec k0,
                            arma::vec a0,
                            arma::vec b0,
                            double mass,
                            double sigma_PY);

void hyper_accelerate_SLI_mv_P(arma::mat mu,
                               arma::mat s2,
                               arma::vec clust,
                               arma::vec &m0,
                               arma::vec &k0,
                               arma::vec a0,
                               arma::vec &b0,
                               arma::vec m1,
                               arma::vec s21,
                               arma::vec tau1,
                               arma::vec tau2,
                               arma::vec a1,
                               arma::vec b1);

void grow_param_SLI_PY_mv_P(arma::mat &mu,
                            arma::mat &s2,
                            arma::vec &v,
                            arma::vec &w,
                            arma::vec u,
                            arma::vec m0,
                            arma::vec k0,
                            arma::vec a0,
                            arma::vec b0,
                            double mass,
                            int n,
                            double sigma_PY);

void grow_param_indep_SLI_PY_mv_P(arma::mat &mu,
                                  arma::mat &s2,
                                  arma::vec &v,
                                  arma::vec &w,
                                  arma::vec &xi,
                                  arma::vec u,
                                  arma::vec m0,
                                  arma::vec k0,
                                  arma::vec a0,
                                  arma::vec b0,
                                  double mass,
                                  int n,
                                  double sigma_PY,
                                  double param_seq_one,
                                  double param_seq_two);

void update_cluster_SLI_mv_P(arma::mat data,
                             arma::mat mu,
                             arma::mat s2,
                             arma::vec &clust,
                             arma::vec w,
                             arma::vec u);

void update_cluster_indep_SLI_mv_P(arma::mat data,
                                   arma::mat mu,
                                   arma::mat s2,
                                   arma::vec &clust,
                                   arma::vec w,
                                   arma::vec xi,
                                   arma::vec u);

/*----------------------------------------------------------------------
 *
 * MIXTURE OF REGRESSION KERNELS
 * LOCATION-SCALE KERNEL
 * SLI functions
 *
 *----------------------------------------------------------------------
 */

void accelerate_SLI_mv_MRK(arma::vec y,
                           arma::mat covs,
                           arma::mat &beta,
                           arma::vec &sigma2,
                           arma::vec &v,
                           arma::vec &w,
                           arma::vec clust,
                           arma::vec beta0,
                           arma::mat Sb0,
                           double a0,
                           double b0,
                           double mass,
                           double sigma_PY);

void hyper_accelerate_SLI_mv_MRK(arma::vec y,
                                 arma::mat covs,
                                 arma::vec clust,
                                 arma::mat beta,
                                 arma::vec sigma2,
                                 arma::vec &beta0,
                                 arma::mat &Sb0,
                                 double a0,
                                 double &b0,
                                 arma::vec beta1,
                                 double k1,
                                 double sb1,
                                 arma::mat Sb1,
                                 double tau1,
                                 double tau2);

void grow_param_SLI_PY_mv_MRK(arma::mat &beta,
                              arma::vec &sigma2,
                              arma::vec &v,
                              arma::vec &w,
                              arma::vec u,
                              arma::vec beta0,
                              arma::mat Sb0,
                              double a0,
                              double b0,
                              double mass,
                              int n,
                              double sigma_PY);

void grow_param_indep_SLI_PY_mv_MRK(arma::mat &beta,
                                    arma::vec &sigma2,
                                    arma::vec &v,
                                    arma::vec &w,
                                    arma::vec &xi,
                                    arma::vec u,
                                    arma::vec beta0,
                                    arma::mat Sb0,
                                    double a0,
                                    double b0,
                                    double mass,
                                    int n,
                                    double sigma_PY,
                                    double param_seq_one,
                                    double param_seq_two);

void update_cluster_SLI_mv_MRK(arma::vec y,
                               arma::mat covs,
                               arma::mat beta,
                               arma::vec sigma2,
                               arma::vec &clust,
                               arma::vec w,
                               arma::vec u);

void update_cluster_indep_SLI_mv_MRK(arma::vec y,
                                     arma::mat covs,
                                     arma::mat beta,
                                     arma::vec sigma2,
                                     arma::vec &clust,
                                     arma::vec w,
                                     arma::vec xi,
                                     arma::vec u);

/*----------------------------------------------------------------------
 *
 * MIXTURE OF REGRESSION KERNELS
 * LOCATION-SCALE KERNEL
 * SLI functions
 *
 *----------------------------------------------------------------------
 */

void accelerate_SLI_mv_MRK_L(arma::vec y,
                             arma::mat covs,
                             arma::mat &beta,
                             double &sigma2,
                             arma::vec &v,
                             arma::vec &w,
                             arma::vec clust,
                             arma::vec beta0,
                             arma::mat Sb0,
                             double a0,
                             double b0,
                             double mass,
                             double sigma_PY);

void hyper_accelerate_SLI_mv_MRK_L(arma::vec y,
                                   arma::mat covs,
                                   arma::vec clust,
                                   arma::mat beta,
                                   arma::vec &beta0,
                                   arma::mat &Sb0,
                                   arma::vec beta1,
                                   double k1,
                                   double sb1,
                                   arma::mat Sb1);

void grow_param_SLI_PY_mv_MRK_L(arma::mat &beta,
                                arma::vec &v,
                                arma::vec &w,
                                arma::vec u,
                                arma::vec beta0,
                                arma::mat Sb0,
                                double mass,
                                int n,
                                double sigma_PY);

void grow_param_indep_SLI_PY_mv_MRK_L(arma::mat &beta,
                                      arma::vec &v,
                                      arma::vec &w,
                                      arma::vec &xi,
                                      arma::vec u,
                                      arma::vec beta0,
                                      arma::mat Sb0,
                                      double mass,
                                      int n,
                                      double sigma_PY,
                                      double param_seq_one,
                                      double param_seq_two);

void update_cluster_SLI_mv_MRK_L(arma::vec y,
                                 arma::mat covs,
                                 arma::mat beta,
                                 double sigma2,
                                 arma::vec &clust,
                                 arma::vec w,
                                 arma::vec u);

void update_cluster_indep_SLI_mv_MRK_L(arma::vec y,
                                       arma::mat covs,
                                       arma::mat beta,
                                       double sigma2,
                                       arma::vec &clust,
                                       arma::vec w,
                                       arma::vec xi,
                                       arma::vec u);

#endif
