// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// optimalAssignment
NumericVector optimalAssignment(NumericMatrix rankings, int leaders, int minGroupSize, int maxGroupSize);
RcppExport SEXP _GroupAssignment_optimalAssignment(SEXP rankingsSEXP, SEXP leadersSEXP, SEXP minGroupSizeSEXP, SEXP maxGroupSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type rankings(rankingsSEXP);
    Rcpp::traits::input_parameter< int >::type leaders(leadersSEXP);
    Rcpp::traits::input_parameter< int >::type minGroupSize(minGroupSizeSEXP);
    Rcpp::traits::input_parameter< int >::type maxGroupSize(maxGroupSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(optimalAssignment(rankings, leaders, minGroupSize, maxGroupSize));
    return rcpp_result_gen;
END_RCPP
}
