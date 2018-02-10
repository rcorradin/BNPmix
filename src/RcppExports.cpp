// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// marginal_DP_multi
Rcpp::List marginal_DP_multi(int nsim, int nburn, int napprox, int nparam, int d, int grid_l, arma::mat data, arma::mat grid, arma::vec conf_start, arma::vec mu_start, arma::mat Lambda_start, double theta, arma::vec m0, arma::mat B0, double nu0, arma::mat sigma, int b1, arma::mat B1, arma::vec m1, arma::mat M1, int s1, arma::mat S1, double t1, double t2, int nupd, int plim, bool FIX);
RcppExport SEXP _BNPmix_marginal_DP_multi(SEXP nsimSEXP, SEXP nburnSEXP, SEXP napproxSEXP, SEXP nparamSEXP, SEXP dSEXP, SEXP grid_lSEXP, SEXP dataSEXP, SEXP gridSEXP, SEXP conf_startSEXP, SEXP mu_startSEXP, SEXP Lambda_startSEXP, SEXP thetaSEXP, SEXP m0SEXP, SEXP B0SEXP, SEXP nu0SEXP, SEXP sigmaSEXP, SEXP b1SEXP, SEXP B1SEXP, SEXP m1SEXP, SEXP M1SEXP, SEXP s1SEXP, SEXP S1SEXP, SEXP t1SEXP, SEXP t2SEXP, SEXP nupdSEXP, SEXP plimSEXP, SEXP FIXSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< int >::type nburn(nburnSEXP);
    Rcpp::traits::input_parameter< int >::type napprox(napproxSEXP);
    Rcpp::traits::input_parameter< int >::type nparam(nparamSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type grid_l(grid_lSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type grid(gridSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type conf_start(conf_startSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_start(mu_startSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Lambda_start(Lambda_startSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type m0(m0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B0(B0SEXP);
    Rcpp::traits::input_parameter< double >::type nu0(nu0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type b1(b1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B1(B1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type m1(m1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type M1(M1SEXP);
    Rcpp::traits::input_parameter< int >::type s1(s1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S1(S1SEXP);
    Rcpp::traits::input_parameter< double >::type t1(t1SEXP);
    Rcpp::traits::input_parameter< double >::type t2(t2SEXP);
    Rcpp::traits::input_parameter< int >::type nupd(nupdSEXP);
    Rcpp::traits::input_parameter< int >::type plim(plimSEXP);
    Rcpp::traits::input_parameter< bool >::type FIX(FIXSEXP);
    rcpp_result_gen = Rcpp::wrap(marginal_DP_multi(nsim, nburn, napprox, nparam, d, grid_l, data, grid, conf_start, mu_start, Lambda_start, theta, m0, B0, nu0, sigma, b1, B1, m1, M1, s1, S1, t1, t2, nupd, plim, FIX));
    return rcpp_result_gen;
END_RCPP
}
// update_cluster_cpp
void update_cluster_cpp(arma::mat data, arma::cube& Lambda, arma::mat& mu, arma::vec& clust, arma::vec& useful, arma::vec m0, arma::mat B0, double nu0, arma::mat sigma, double theta, int napprox);
RcppExport SEXP _BNPmix_update_cluster_cpp(SEXP dataSEXP, SEXP LambdaSEXP, SEXP muSEXP, SEXP clustSEXP, SEXP usefulSEXP, SEXP m0SEXP, SEXP B0SEXP, SEXP nu0SEXP, SEXP sigmaSEXP, SEXP thetaSEXP, SEXP napproxSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type Lambda(LambdaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type clust(clustSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type useful(usefulSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type m0(m0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B0(B0SEXP);
    Rcpp::traits::input_parameter< double >::type nu0(nu0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type napprox(napproxSEXP);
    update_cluster_cpp(data, Lambda, mu, clust, useful, m0, B0, nu0, sigma, theta, napprox);
    return R_NilValue;
END_RCPP
}
// update_hyperparameters
double update_hyperparameters(int n, double& theta, arma::cube Lambda, arma::mat mu, arma::vec clust, arma::vec useful, arma::vec& m0, arma::mat& B0, double nu0, arma::mat& sigma, int b1, arma::mat B1, arma::vec m1, arma::mat M1, int s1, arma::mat S1, double t1, double t2, bool FIX);
RcppExport SEXP _BNPmix_update_hyperparameters(SEXP nSEXP, SEXP thetaSEXP, SEXP LambdaSEXP, SEXP muSEXP, SEXP clustSEXP, SEXP usefulSEXP, SEXP m0SEXP, SEXP B0SEXP, SEXP nu0SEXP, SEXP sigmaSEXP, SEXP b1SEXP, SEXP B1SEXP, SEXP m1SEXP, SEXP M1SEXP, SEXP s1SEXP, SEXP S1SEXP, SEXP t1SEXP, SEXP t2SEXP, SEXP FIXSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Lambda(LambdaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type clust(clustSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type useful(usefulSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type m0(m0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type B0(B0SEXP);
    Rcpp::traits::input_parameter< double >::type nu0(nu0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type b1(b1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B1(B1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type m1(m1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type M1(M1SEXP);
    Rcpp::traits::input_parameter< int >::type s1(s1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S1(S1SEXP);
    Rcpp::traits::input_parameter< double >::type t1(t1SEXP);
    Rcpp::traits::input_parameter< double >::type t2(t2SEXP);
    Rcpp::traits::input_parameter< bool >::type FIX(FIXSEXP);
    rcpp_result_gen = Rcpp::wrap(update_hyperparameters(n, theta, Lambda, mu, clust, useful, m0, B0, nu0, sigma, b1, B1, m1, M1, s1, S1, t1, t2, FIX));
    return rcpp_result_gen;
END_RCPP
}
// update_parameters
void update_parameters(arma::mat data, arma::cube& Lambda, arma::mat& mu, arma::vec& clust, arma::vec useful, arma::vec m0, arma::mat B0, double nu0, arma::mat sigma);
RcppExport SEXP _BNPmix_update_parameters(SEXP dataSEXP, SEXP LambdaSEXP, SEXP muSEXP, SEXP clustSEXP, SEXP usefulSEXP, SEXP m0SEXP, SEXP B0SEXP, SEXP nu0SEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type Lambda(LambdaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type clust(clustSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type useful(usefulSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type m0(m0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B0(B0SEXP);
    Rcpp::traits::input_parameter< double >::type nu0(nu0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    update_parameters(data, Lambda, mu, clust, useful, m0, B0, nu0, sigma);
    return R_NilValue;
END_RCPP
}
// find_part
Rcpp::List find_part(arma::mat clust_mat);
RcppExport SEXP _BNPmix_find_part(SEXP clust_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type clust_mat(clust_matSEXP);
    rcpp_result_gen = Rcpp::wrap(find_part(clust_mat));
    return rcpp_result_gen;
END_RCPP
}
// pairwise_mat
Rcpp::List pairwise_mat(arma::mat clust_mat);
RcppExport SEXP _BNPmix_pairwise_mat(SEXP clust_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type clust_mat(clust_matSEXP);
    rcpp_result_gen = Rcpp::wrap(pairwise_mat(clust_mat));
    return rcpp_result_gen;
END_RCPP
}
// est_ISE_2D
double est_ISE_2D(arma::vec estimated, arma::vec teoric, arma::mat grid);
RcppExport SEXP _BNPmix_est_ISE_2D(SEXP estimatedSEXP, SEXP teoricSEXP, SEXP gridSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type estimated(estimatedSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type teoric(teoricSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type grid(gridSEXP);
    rcpp_result_gen = Rcpp::wrap(est_ISE_2D(estimated, teoric, grid));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BNPmix_marginal_DP_multi", (DL_FUNC) &_BNPmix_marginal_DP_multi, 27},
    {"_BNPmix_update_cluster_cpp", (DL_FUNC) &_BNPmix_update_cluster_cpp, 11},
    {"_BNPmix_update_hyperparameters", (DL_FUNC) &_BNPmix_update_hyperparameters, 19},
    {"_BNPmix_update_parameters", (DL_FUNC) &_BNPmix_update_parameters, 9},
    {"_BNPmix_find_part", (DL_FUNC) &_BNPmix_find_part, 1},
    {"_BNPmix_pairwise_mat", (DL_FUNC) &_BNPmix_pairwise_mat, 1},
    {"_BNPmix_est_ISE_2D", (DL_FUNC) &_BNPmix_est_ISE_2D, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_BNPmix(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
