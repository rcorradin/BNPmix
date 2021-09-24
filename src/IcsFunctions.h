#ifndef cICS_H
#define cICS_H

/*----------------------------------------------------------------------
 *
 * UNIVARIATE
 * LOCATION
 * ICS functions
 *
 *----------------------------------------------------------------------
 */

void accelerate_ICS_L(arma::vec data,
                      arma::vec &mu,
                      double &s2,
                      arma::vec clust,
                      double m0,
                      double s20,
                      double a0,
                      double b0);

void hyper_accelerate_ICS_L(arma::vec mu,
                            double &m0,
                            double &s20,
                            double m1,
                            double k1,
                            double a1,
                            double b1);

void para_clean_ICS_L(arma::vec &mu,
                      arma::vec &clust);

void para_clean_ICS_L_export(arma::vec &mu,
                             arma::vec &mujoin,
                             arma::vec &probjoin,
                             arma::vec &clust);

void simu_trunc_PY_L(arma::vec &mutemp,
                     arma::vec &freqtemp,
                     double mass,
                     double m0,
                     double s20,
                     int napprox,
                     double sigma_PY);

void clust_update_ICS_L(arma::vec data,
                        arma::vec mujoin,
                        double s2,
                        arma::vec probjoin,
                        arma::vec &clust);

/*----------------------------------------------------------------------
 *
 * UNIVARIATE
 * LOCATION SCALE
 * ICS functions
 *
 *----------------------------------------------------------------------
 */

void accelerate_ICS(arma::vec data,
                    arma::vec &mu,
                    arma::vec &s2,
                    arma::vec clust,
                    double m0,
                    double k0,
                    double a0,
                    double b0);

void hyper_accelerate_ICS(arma::vec mu,
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

void para_clean_ICS(arma::vec &mu,
                    arma::vec &s2,
                    arma::vec &clust);

void para_clean_ICS_export(arma::vec &mu,
                           arma::vec &s2,
                           arma::vec &mujoin,
                           arma::vec &s2join,
                           arma::vec &probjoin,
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

void clust_update_ICS(arma::vec data,
                      arma::vec mujoin,
                      arma::vec s2join,
                      arma::vec probjoin,
                      arma::vec &clust);

/*----------------------------------------------------------------------
 *
 * MULTIVARIATE
 * LOCATION
 * ICS functions
 *
 *----------------------------------------------------------------------
 */

void accelerate_ICS_mv_L(arma::mat data,
                         arma::mat &mu,
                         arma::mat &s2,
                         arma::vec clust,
                         arma::vec m0,
                         arma::mat S20,
                         arma::mat S0,
                         double n0);

void hyper_accelerate_ICS_mv_L(arma::mat mu,
                               arma::vec &m0,
                               arma::mat &S20,
                               arma::vec m1,
                               double k1,
                               double theta1,
                               arma::mat Theta1);

void para_clean_ICS_mv_L(arma::mat &mu,
                         arma::vec &clust);

void para_clean_ICS_mv_L_export(arma::mat &mu,
                                arma::mat &mujoin,
                                arma::vec &probjoin,
                                arma::vec &clust);

void simu_trunc_PY_mv_L(arma::mat &mutemp,
                        arma::vec &freqtemp,
                        double mass,
                        arma::vec m0,
                        arma::mat S20,
                        int napprox,
                        double sigma_PY);

void clust_update_ICS_mv_L(arma::mat data,
                           arma::mat mujoin,
                           arma::mat s2,
                           arma::vec probjoin,
                           arma::vec &clust);

/*----------------------------------------------------------------------
 *
 * MULTIVARIATE
 * LOCATION SCALE
 * ICS functions
 *
 *----------------------------------------------------------------------
 */

void accelerate_ICS_mv(arma::mat data,
                       arma::mat &mu,
                       arma::cube &s2,
                       arma::vec clust,
                       arma::vec m0,
                       double k0,
                       arma::mat S0,
                       double n0);

void hyper_accelerate_ICS_mv_LS(arma::mat mu,
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

void para_clean_ICS_mv(arma::mat &mu,
                       arma::cube &s2,
                       arma::vec &clust);

void para_clean_ICS_mv_export(arma::mat &mu,
                              arma::cube &s2,
                              arma::mat &mujoin,
                              arma::cube &s2join,
                              arma::vec &probjoin,
                              arma::vec &clust);

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

void clust_update_ICS_mv(arma::mat data,
                         arma::mat mujoin,
                         arma::cube s2join,
                         arma::vec probjoin,
                         arma::vec &clust,
                         //-*-temp
                         double &new_clust);

/*----------------------------------------------------------------------
 *
 * MULTIVARIATE
 * PRODUCT KERNEL
 * ICS functions
 *
 *----------------------------------------------------------------------
 */

void accelerate_ICS_mv_P(arma::mat data,
                         arma::mat &mu,
                         arma::mat &s2,
                         arma::vec clust,
                         arma::vec m0,
                         arma::vec k0,
                         arma::vec a0,
                         arma::vec b0);

void hyper_accelerate_ICS_mv_P(arma::mat mu,
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

void para_clean_ICS_mv_P(arma::mat &mu,
                         arma::mat &s2,
                         arma::vec &clust);

void para_clean_ICS_mv_P_export(arma::mat &mu,
                                arma::mat &s2,
                                arma::mat &mujoin,
                                arma::mat &s2join,
                                arma::vec &probjoin,
                                arma::vec &clust);

void simu_trunc_PY_mv_P(arma::mat &mutemp,
                        arma::mat &s2temp,
                        arma::vec &freqtemp,
                        double mass,
                        arma::vec m0,
                        arma::vec k0,
                        arma::vec a0,
                        arma::vec b0,
                        int napprox,
                        double sigma_PY);

void clust_update_ICS_mv_P(arma::mat data,
                           arma::mat mujoin,
                           arma::mat s2join,
                           arma::vec probjoin,
                           arma::vec &clust);

/*----------------------------------------------------------------------
 *
 * MIXTURE OF REGRESSION KERNELS
 * LOCATION-SCALE KERNEL
 * ICS functions
 *
 *----------------------------------------------------------------------
 */

void accelerate_ICS_mv_MRK(arma::vec y,
                           arma::mat covs,
                           arma::mat &beta,
                           arma::vec &sigma2,
                           arma::vec clust,
                           arma::vec beta0,
                           arma::mat Sb0,
                           double a0,
                           double b0);

void hyper_accelerate_ICS_mv_MRK(arma::vec y,
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

void para_clean_ICS_mv_MRK(arma::mat &beta,
                           arma::vec &sigma2,
                           arma::vec &clust);

void para_clean_ICS_mv_MRK_export(arma::mat &beta,
                                  arma::vec &sigma2,
                                  arma::mat &betajoin,
                                  arma::vec &sigma2join,
                                  arma::vec &probjoin,
                                  arma::vec &clust);

void simu_trunc_PY_mv_MRK(arma::mat &betatemp,
                          arma::vec &sigma2temp,
                          arma::vec &freqtemp,
                          double mass,
                          arma::vec beta0,
                          arma::mat Sb0,
                          double a0,
                          double b0,
                          int napprox,
                          double sigma_PY);

void clust_update_ICS_mv_MRK(arma::vec y,
                             arma::mat covs,
                             arma::mat betajoin,
                             arma::vec sigma2join,
                             arma::vec probjoin,
                             arma::vec &clust);

/*----------------------------------------------------------------------
 *
 * MIXTURE OF REGRESSION KERNELS - LOCATION
 * LOCATION-SCALE KERNEL
 * ICS functions
 *
 *----------------------------------------------------------------------
 */

void accelerate_ICS_mv_MRK_L(arma::vec y,
                             arma::mat covs,
                             arma::mat &beta,
                             double &sigma2,
                             arma::vec clust,
                             arma::vec beta0,
                             arma::mat Sb0,
                             double a0,
                             double b0);

void hyper_accelerate_ICS_mv_MRK_L(arma::vec y,
                                   arma::mat covs,
                                   arma::vec clust,
                                   arma::mat beta,
                                   arma::vec &beta0,
                                   arma::mat &Sb0,
                                   double a0,
                                   double &b0,
                                   arma::vec beta1,
                                   double k1,
                                   double sb1,
                                   arma::mat Sb1);

void para_clean_ICS_mv_MRK_L(arma::mat &beta,
                             arma::vec &clust);

void para_clean_ICS_mv_MRK_L_export(arma::mat &beta,
                                    arma::mat &betajoin,
                                    arma::vec &probjoin,
                                    arma::vec &clust);

void simu_trunc_PY_mv_MRK_L(arma::mat &betatemp,
                            arma::vec &freqtemp,
                            double mass,
                            arma::vec beta0,
                            arma::mat Sb0,
                            int napprox,
                            double sigma_PY);

void clust_update_ICS_mv_MRK_L(arma::vec y,
                               arma::mat covs,
                               arma::mat betajoin,
                               double sigma2,
                               arma::vec probjoin,
                               arma::vec &clust);

#endif
