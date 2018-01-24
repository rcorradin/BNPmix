#ifndef UPDATE_CLUSTER_H
#define UPDATE_CLUSTER_H
void update_cluster_cpp(arma::mat data, 
                        arma::cube& Lambda, 
                        arma::mat& mu, 
                        arma::vec& clust, 
                        arma::vec &useful, 
                        arma::vec m0, 
                        arma::mat B0, 
                        double nu0, 
                        arma::mat sigma, 
                        double theta, 
                        int napprox);
#endif