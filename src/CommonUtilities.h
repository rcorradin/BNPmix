#ifndef COMMON_H
#define COMMON_H

arma::vec freq_vec(arma::vec vector);

arma::vec eval_density_L(arma::vec grid,
                         arma::vec mu,
                         double s2,
                         arma::vec probs);

arma::vec eval_density(arma::vec grid,
                       arma::vec mu,
                       arma::vec s2,
                       arma::vec probs);

arma::vec eval_density_MAR(arma::vec grid,
                           arma::vec mu,
                           arma::vec s2,
                           arma::vec probs,
                           double m0,
                           double k0,
                           double a0,
                           double b0);

arma::vec eval_density_mv(arma::mat grid,
                          arma::mat mu,
                          arma::cube s2,
                          arma::vec probs);

arma::vec eval_density_mv_L(arma::mat grid,
                            arma::mat mu,
                            arma::mat s2,
                            arma::vec probs);

arma::vec eval_density_mv_P(arma::mat grid,
                            arma::mat mu,
                            arma::mat s2,
                            arma::vec probs);

arma::mat eval_density_mv_MKR(arma::mat grid_covs,
                              arma::vec grid_response,
                              arma::mat beta,
                              arma::vec sigma2,
                              arma::vec prob);

arma::mat eval_density_mv_MKR_L(arma::mat grid_covs,
                                arma::vec grid_response,
                                arma::mat beta,
                                double sigma2,
                                arma::vec prob);
#endif
