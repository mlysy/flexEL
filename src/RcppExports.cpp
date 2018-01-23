// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// lambdaNR
Rcpp::List lambdaNR(Eigen::MatrixXd G, Eigen::VectorXd lambda0, int maxIter, double eps, bool verbose);
RcppExport SEXP _bayesEL_lambdaNR(SEXP GSEXP, SEXP lambda0SEXP, SEXP maxIterSEXP, SEXP epsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type G(GSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(lambdaNR(G, lambda0, maxIter, eps, verbose));
    return rcpp_result_gen;
END_RCPP
}
// logELmean
double logELmean(int nObs, int nEqs, Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::VectorXd beta, Eigen::VectorXd lambda0, int maxIter, double eps);
RcppExport SEXP _bayesEL_logELmean(SEXP nObsSEXP, SEXP nEqsSEXP, SEXP ySEXP, SEXP XSEXP, SEXP betaSEXP, SEXP lambda0SEXP, SEXP maxIterSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nObs(nObsSEXP);
    Rcpp::traits::input_parameter< int >::type nEqs(nEqsSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(logELmean(nObs, nEqs, y, X, beta, lambda0, maxIter, eps));
    return rcpp_result_gen;
END_RCPP
}
// MeanReg_post
Eigen::MatrixXd MeanReg_post(int nObs, int nEqs, Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::VectorXd lambda0, int nsamples, int nburn, Eigen::VectorXd betaInit, Eigen::VectorXd sigs, int maxIter, double eps);
RcppExport SEXP _bayesEL_MeanReg_post(SEXP nObsSEXP, SEXP nEqsSEXP, SEXP ySEXP, SEXP XSEXP, SEXP lambda0SEXP, SEXP nsamplesSEXP, SEXP nburnSEXP, SEXP betaInitSEXP, SEXP sigsSEXP, SEXP maxIterSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nObs(nObsSEXP);
    Rcpp::traits::input_parameter< int >::type nEqs(nEqsSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< int >::type nsamples(nsamplesSEXP);
    Rcpp::traits::input_parameter< int >::type nburn(nburnSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type betaInit(betaInitSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type sigs(sigsSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(MeanReg_post(nObs, nEqs, y, X, lambda0, nsamples, nburn, betaInit, sigs, maxIter, eps));
    return rcpp_result_gen;
END_RCPP
}
// logELquant
double logELquant(int nObs, int nEqs, Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::VectorXd beta, double alpha, Eigen::VectorXd lambda0, int maxIter, double eps);
RcppExport SEXP _bayesEL_logELquant(SEXP nObsSEXP, SEXP nEqsSEXP, SEXP ySEXP, SEXP XSEXP, SEXP betaSEXP, SEXP alphaSEXP, SEXP lambda0SEXP, SEXP maxIterSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nObs(nObsSEXP);
    Rcpp::traits::input_parameter< int >::type nEqs(nEqsSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(logELquant(nObs, nEqs, y, X, beta, alpha, lambda0, maxIter, eps));
    return rcpp_result_gen;
END_RCPP
}
// QuantReg_post
Eigen::MatrixXd QuantReg_post(int nObs, int nEqs, Eigen::VectorXd y, Eigen::MatrixXd X, double alpha, Eigen::VectorXd lambda0, int nsamples, int nburn, Eigen::VectorXd betaInit, Eigen::VectorXd sigs, int maxIter, double eps);
RcppExport SEXP _bayesEL_QuantReg_post(SEXP nObsSEXP, SEXP nEqsSEXP, SEXP ySEXP, SEXP XSEXP, SEXP alphaSEXP, SEXP lambda0SEXP, SEXP nsamplesSEXP, SEXP nburnSEXP, SEXP betaInitSEXP, SEXP sigsSEXP, SEXP maxIterSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nObs(nObsSEXP);
    Rcpp::traits::input_parameter< int >::type nEqs(nEqsSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< int >::type nsamples(nsamplesSEXP);
    Rcpp::traits::input_parameter< int >::type nburn(nburnSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type betaInit(betaInitSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type sigs(sigsSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(QuantReg_post(nObs, nEqs, y, X, alpha, lambda0, nsamples, nburn, betaInit, sigs, maxIter, eps));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bayesEL_lambdaNR", (DL_FUNC) &_bayesEL_lambdaNR, 5},
    {"_bayesEL_logELmean", (DL_FUNC) &_bayesEL_logELmean, 8},
    {"_bayesEL_MeanReg_post", (DL_FUNC) &_bayesEL_MeanReg_post, 11},
    {"_bayesEL_logELquant", (DL_FUNC) &_bayesEL_logELquant, 9},
    {"_bayesEL_QuantReg_post", (DL_FUNC) &_bayesEL_QuantReg_post, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_bayesEL(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
