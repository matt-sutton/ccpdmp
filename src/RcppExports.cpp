// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// linear_inv_t
arma::vec linear_inv_t(double a, double b, double u, double tmax);
RcppExport SEXP _ccpdmp_linear_inv_t(SEXP aSEXP, SEXP bSEXP, SEXP uSEXP, SEXP tmaxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type u(uSEXP);
    Rcpp::traits::input_parameter< double >::type tmax(tmaxSEXP);
    rcpp_result_gen = Rcpp::wrap(linear_inv_t(a, b, u, tmax));
    return rcpp_result_gen;
END_RCPP
}
// exp_inv_t
double exp_inv_t(double a, double b, double u);
RcppExport SEXP _ccpdmp_exp_inv_t(SEXP aSEXP, SEXP bSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(exp_inv_t(a, b, u));
    return rcpp_result_gen;
END_RCPP
}
// sim_rate_poly
List sim_rate_poly(arma::vec eval_times, arma::vec eval_rates, int poly_order);
RcppExport SEXP _ccpdmp_sim_rate_poly(SEXP eval_timesSEXP, SEXP eval_ratesSEXP, SEXP poly_orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type eval_times(eval_timesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type eval_rates(eval_ratesSEXP);
    Rcpp::traits::input_parameter< int >::type poly_order(poly_orderSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_rate_poly(eval_times, eval_rates, poly_order));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ccpdmp_linear_inv_t", (DL_FUNC) &_ccpdmp_linear_inv_t, 4},
    {"_ccpdmp_exp_inv_t", (DL_FUNC) &_ccpdmp_exp_inv_t, 3},
    {"_ccpdmp_sim_rate_poly", (DL_FUNC) &_ccpdmp_sim_rate_poly, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_ccpdmp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
