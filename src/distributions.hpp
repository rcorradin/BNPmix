#ifndef GAUSSIAN_H
#define GAUSSIAN_H

//  / ___| __ _ _   _ ___ ___(_) __ _ _ __
// | |  _ / _` | | | / __/ __| |/ _` | '_ \
// | |_| | (_| | |_| \__ \__ \ | (_| | | | |
//  \____|\__,_|\__,_|___/___/_|\__,_|_| |_|

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

//  _       ____  _             _            _
// | |_    / ___|| |_ _   _  __| | ___ _ __ | |_
// | __|___\___ \| __| | | |/ _` |/ _ \ '_ \| __|
// | ||_____|__) | |_| |_| | (_| |  __/ | | | |_
//  \__|   |____/ \__|\__,_|\__,_|\___|_| |_|\__|

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

// __        ___     _                _
// \ \      / (_)___| |__   __ _ _ __| |_
//  \ \ /\ / /| / __| '_ \ / _` | '__| __|
//   \ V  V / | \__ \ | | | (_| | |  | |_
//    \_/\_/  |_|___/_| |_|\__,_|_|   \__|

arma::mat rWishartMat(int df,
                      arma::mat lambda);

//  ____  _                   _
// |  _ \(_)___  ___ _ __ ___| |_ ___
// | | | | / __|/ __| '__/ _ \ __/ _ \
// | |_| | \__ \ (__| | |  __/ ||  __/
// |____/|_|___/\___|_|  \___|\__\___|
// | | | |_ __ (_)/ _| ___  _ __ _ __ ___
// | | | | '_ \| | |_ / _ \| '__| '_ ` _ \
// | |_| | | | | |  _| (_) | |  | | | | | |
//  \___/|_| |_|_|_|  \___/|_|  |_| |_| |_|

int rintnunif(arma::vec prob,
              int a);

#endif
