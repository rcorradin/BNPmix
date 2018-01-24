#ifndef GAUSSIAN_H
#define GAUSSIAN_H

double dmvnrm_ar(arma::rowvec x, 
                 arma::rowvec mean,  
                 arma::mat sigma, 
                 bool logd = false);

double dmvnrm_prec(arma::rowvec x, 
                   arma::rowvec mean,  
                   arma::mat lambda, 
                   bool logd = false);

arma::vec dmvnrm_ar_mat(arma::mat x,  
                          arma::rowvec mean,  
                          arma::mat sigma, 
                          bool logd = false);

arma::vec dmvnrm_prec_mat(arma::mat x,  
                          arma::rowvec mean,  
                          arma::mat lambda, 
                          bool logd = false);

double dmn_prec_det(arma::rowvec x, 
                    arma::rowvec mean, 
                    arma::mat lambda, 
                    double det);

arma::mat rmvnormMat(int n, 
					 arma::vec mu, 
					 arma::mat sigma);
#endif