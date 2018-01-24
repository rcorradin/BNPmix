#ifndef UPDATE_PARAMETERS_H
#define UPDATE_PARAMETERS_H
void update_parameters(arma::mat data, 
                       arma::cube& Lambda, 
                       arma::mat& mu, 
                       arma::vec& clust, 
                       arma::vec useful, 
                       arma::vec m0, 
                       arma::mat B0, 
                       double nu0, 
                       arma::mat sigma);
#endif