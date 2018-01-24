#ifndef UPDATE_HYPERPARAMETERS_H
#define UPDATE_HYPERPARAMETERS_H
double update_hyperparameters(int n, 
                              double& theta, 
                              arma::cube Lambda, 
                              arma::mat mu, 
                              arma::vec clust, 
                              arma::vec useful, 
                              arma::vec& m0, 
                              arma::mat& B0, 
                              double nu0, 
                              arma::mat& sigma, 
                              int b1, 
                              arma::mat B1, 
                              arma::vec m1, 
                              arma::mat M1, 
                              int s1, 
                              arma::mat S1, 
                              double t1, 
                              double t2);
#endif