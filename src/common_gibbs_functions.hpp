/*
  common gibbs function for different models
*/

void para_cleanser(arma::cube &Lambda,
                   arma::mat &mu,
                   arma::vec &clust);

void update_cluster_cpp(arma::mat data,
                        arma::cube &Lambda,
                        arma::mat &mu,
                        arma::vec &clust,
                        arma::vec m0,
                        arma::mat B0,
                        double nu0,
                        arma::mat sigma,
                        double theta,
                        int napprox);

void update_parameters(arma::mat data,
                       arma::cube &Lambda,
                       arma::mat &mu,
                       arma::vec &clust,
                       arma::vec m0,
                       arma::mat B0,
                       double nu0,
                       arma::mat sigma);

arma::vec update_distribution(arma::mat grid,
                              int grid_l,
                              arma::mat mu,
                              arma::cube Lambda,
                              arma::vec clust,
                              double theta);
