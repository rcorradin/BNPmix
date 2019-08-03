#ifndef MAR_H
#define MAR_H

/*----------------------------------------------------------------------
 *
 * UNIVARIATE
 * LOCATION KERNEL
 * ICS functions
 *
 *----------------------------------------------------------------------
 */

void accelerate_MAR_L(arma::vec data,
                      arma::vec &mu,
                      double s2,
                      arma::vec clust,
                      double m0,
                      double s20,
                      double a0,
                      double b0);

void hyper_accelerate_MAR_L(arma::vec mu,
                            double &m0,
                            double &s20,
                            double m1,
                            double k1,
                            double a1,
                            double b1);

void para_clean_MAR_L(arma::vec &mu,
                      arma::vec &clust);

void clust_update_MAR_PY_L(arma::vec data,
                           arma::vec &mu,
                           double &s2,
                           arma::vec &clust,
                           double mass,
                           double m0,
                           double s20,
                           double a0,
                           double b0,
                           double sigma_PY);

/*----------------------------------------------------------------------
 *
 * UNIVARIATE
 * LOCATION-SCALE KERNEL
 * ICS functions
 *
 *----------------------------------------------------------------------
 */

void accelerate_MAR(arma::vec data,
                    arma::vec &mu,
                    arma::vec &s2,
                    arma::vec clust,
                    double m0,
                    double k0,
                    double a0,
                    double b0);

void hyper_accelerate_MAR(arma::vec mu,
                          arma::vec s2,
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

void para_clean_MAR(arma::vec &mu,
                    arma::vec &s2,
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
                         double sigma_PY);

/*----------------------------------------------------------------------
 *
 * MULTIVARIATE
 * LOCATION KERNEL
 * ICS functions
 *
 *----------------------------------------------------------------------
 */

void accelerate_MAR_mv_L(arma::mat data,
                         arma::mat &mu,
                         arma::mat &s2,
                         arma::vec clust,
                         arma::vec m0,
                         arma::mat S20,
                         arma::mat S0,
                         double n0);

void hyper_accelerate_MAR_mv_L(arma::mat mu,
                               arma::vec &m0,
                               arma::mat &S20,
                               arma::vec m1,
                               double k1,
                               double theta1,
                               arma::mat Theta1);

void para_clean_MAR_mv_L(arma::mat &mu,
                         arma::vec &clust);

void clust_update_MAR_PY_mv_L(arma::mat data,
                              arma::mat &mu,
                              arma::mat s2,
                              arma::vec &clust,
                              double mass,
                              arma::vec m0,
                              arma::mat S20,
                              arma::mat S0,
                              double n0,
                              double sigma_PY);

/*----------------------------------------------------------------------
 *
 * MULTIVARIATE
 * LOCATION-SCALE KERNEL
 * ICS functions
 *
 *----------------------------------------------------------------------
 */

void accelerate_MAR_mv(arma::mat data,
                       arma::mat &mu,
                       arma::cube &s2,
                       arma::vec clust,
                       arma::vec m0,
                       double k0,
                       arma::mat S0,
                       double n0);

void hyper_accelerate_MAR_mv_LS(arma::mat mu,
                                arma::cube s2,
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

void para_clean_MAR_mv(arma::mat &mu,
                       arma::cube &s2,
                       arma::vec &clust);

void clust_update_MAR_PY_mv(arma::mat data,
                            arma::mat &mu,
                            arma::cube &s2,
                            arma::vec &clust,
                            double mass,
                            arma::vec m0,
                            double k0,
                            arma::mat S0,
                            double n0,
                            double sigma_PY);

/*----------------------------------------------------------------------
 *
 * MULTIVARIATE
 * PRODUCT KERNEL
 * MAR functions
 *
 *----------------------------------------------------------------------
 */

void accelerate_MAR_mv_P(arma::mat data,
                         arma::mat &mu,
                         arma::mat &s2,
                         arma::vec clust,
                         arma::vec m0,
                         arma::vec k0,
                         arma::vec a0,
                         arma::vec b0);

void hyper_accelerate_MAR_mv_P(arma::mat mu,
                               arma::mat s2,
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

void para_clean_MAR_mv_P(arma::mat &mu,
                         arma::mat &s2,
                         arma::vec &clust);

void clust_update_MAR_PY_mv_P(arma::mat data,
                              arma::mat &mu,
                              arma::mat &s2,
                              arma::vec &clust,
                              double mass,
                              arma::vec m0,
                              arma::vec k0,
                              arma::vec a0,
                              arma::vec b0,
                              double sigma_PY);

/*----------------------------------------------------------------------
 *
 * MIXTURE OF REGRESSION KERNELS
 * LOCATION-SCALE KERNEL
 * MAR functions
 *
 *----------------------------------------------------------------------
 */

void accelerate_MAR_mv_MRK(arma::vec y,
                           arma::mat covs,
                           arma::mat &beta,
                           arma::vec &sigma2,
                           arma::vec clust,
                           arma::vec beta0,
                           arma::mat Sb0,
                           double a0,
                           double b0);

void hyper_accelerate_MAR_mv_MRK(arma::vec y,
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

void para_clean_MAR_mv_MRK(arma::mat &beta,
                           arma::vec &sigma2,
                           arma::vec &clust);

void clust_update_MAR_mv_MRK(arma::vec y,
                             arma::mat covs,
                             arma::mat &beta,
                             arma::vec &sigma2,
                             arma::vec &clust,
                             double mass,
                             arma::vec beta0,
                             arma::mat Sb0,
                             double a0,
                             double b0,
                             double sigma_PY,
                             int napprox);

#endif
