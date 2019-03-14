#ifndef COMMON_H
#define COMMON_H

arma::vec freq_vec(arma::vec vector);

arma::vec eval_density(arma::vec grid,
                       arma::vec mu,
                       arma::vec s2,
                       arma::vec probs);

arma::vec eval_density_mv(arma::mat grid,
                          arma::mat mu,
                          arma::cube s2,
                          arma::vec probs);

arma::vec cond_dist(arma::mat mu,
                    arma::cube s2,
                    arma::vec probs,
                    arma::vec uniquey,
                    double upperbound);

#endif
