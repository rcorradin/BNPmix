#ifndef TSTUDENT_H
#define TSTUDENT_H

double dmvt_ar(arma::rowvec x, 
               arma::rowvec mean,  
               arma::mat sigma, 
               int df, 
               bool logd = false);

double dmvt_prec(arma::rowvec x, 
                 arma::rowvec mean,  
                 arma::mat lambda, 
                 int df, 
                 bool logd = false);

arma::vec dmvt_prec_mat(arma::mat x, 
                 arma::rowvec mean,  
                 arma::mat lambda, 
                 int df, 
                 bool logd = false);


#endif