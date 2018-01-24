#ifndef UPDATE_DISTRIBUTION_H
#define UPDATE_DISTRIBUTION_H
arma::vec update_distribution(arma::mat grid, 
                              int grid_l,
                              arma::mat mu, 
                              arma::cube Lambda, 
                              arma::vec useful, 
                              arma::vec clust, 
                              double theta, 
                              int n);
#endif