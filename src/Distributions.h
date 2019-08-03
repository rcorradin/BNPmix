#ifndef GAUSSIAN_H
#define GAUSSIAN_H

int rintnunif(arma::vec weights);

int rintnunif_log(arma::vec lweights);

int rintnunifw(arma::vec freq,
               double mass);

arma::vec rdirich_mass(arma::vec freq,
                       double mass);

arma::vec rdirich_mass_tot(arma::vec freq,
                 double mass);

double dt_ls(double x,
             double gdl,
             double mu,
             double sigma);

double dt_ls_mv(arma::vec x,
                double df,
                arma::vec mean,
                arma::mat sigma);

#endif
