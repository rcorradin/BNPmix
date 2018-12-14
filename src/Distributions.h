#ifndef GAUSSIAN_H
#define GAUSSIAN_H

// NON-Uniform discrete distribution
int rintnunif(arma::vec weights);

// NON-Uniform discrete distribution - log
int rintnunif_log(arma::vec lweights);

// normalized freq + mass
int rintnunifw(arma::vec freq,
               double mass);

// Dirichlet
arma::vec rdirich_mass(arma::vec freq,
                       double mass);

// Dirichlet + sum gamma
arma::vec rdirich_mass_tot(arma::vec freq,
                 double mass);

// univariate t student density
double dt_ls(double x,
             double gdl,
             double mu,
             double sigma);

// multivariate t student density
double dt_ls_mv(arma::vec x,
                double df,
                arma::vec mean,
                arma::mat sigma);

#endif
