// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// lambdaNR
Rcpp::List lambdaNR(Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::MatrixXd G, Eigen::VectorXd lambda0, int nObs, int nEqs, int maxIter, double eps, bool verbose);
RcppExport SEXP _bayesELnew_lambdaNR(SEXP ySEXP, SEXP XSEXP, SEXP GSEXP, SEXP lambda0SEXP, SEXP nObsSEXP, SEXP nEqsSEXP, SEXP maxIterSEXP, SEXP epsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type G(GSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< int >::type nObs(nObsSEXP);
    Rcpp::traits::input_parameter< int >::type nEqs(nEqsSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(lambdaNR(y, X, G, lambda0, nObs, nEqs, maxIter, eps, verbose));
    return rcpp_result_gen;
END_RCPP
}
// logEL_MeanReg
double logEL_MeanReg(Eigen::VectorXd y, Eigen::MatrixXd X, int nObs, int nEqs, Eigen::VectorXd beta, Eigen::VectorXd lambda0, int maxIter, double eps);
RcppExport SEXP _bayesELnew_logEL_MeanReg(SEXP ySEXP, SEXP XSEXP, SEXP nObsSEXP, SEXP nEqsSEXP, SEXP betaSEXP, SEXP lambda0SEXP, SEXP maxIterSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type nObs(nObsSEXP);
    Rcpp::traits::input_parameter< int >::type nEqs(nEqsSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(logEL_MeanReg(y, X, nObs, nEqs, beta, lambda0, maxIter, eps));
    return rcpp_result_gen;
END_RCPP
}
// PostSample_MeanReg
Eigen::MatrixXd PostSample_MeanReg(Eigen::VectorXd y, Eigen::MatrixXd X, int nObs, int nEqs, Eigen::VectorXd lambda0, int nsamples, int nburn, Eigen::VectorXd betaInit, Eigen::VectorXd sigs, int maxIter, double eps);
RcppExport SEXP _bayesELnew_PostSample_MeanReg(SEXP ySEXP, SEXP XSEXP, SEXP nObsSEXP, SEXP nEqsSEXP, SEXP lambda0SEXP, SEXP nsamplesSEXP, SEXP nburnSEXP, SEXP betaInitSEXP, SEXP sigsSEXP, SEXP maxIterSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type nObs(nObsSEXP);
    Rcpp::traits::input_parameter< int >::type nEqs(nEqsSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< int >::type nsamples(nsamplesSEXP);
    Rcpp::traits::input_parameter< int >::type nburn(nburnSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type betaInit(betaInitSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type sigs(sigsSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(PostSample_MeanReg(y, X, nObs, nEqs, lambda0, nsamples, nburn, betaInit, sigs, maxIter, eps));
    return rcpp_result_gen;
END_RCPP
}
// logEL_QuantReg
double logEL_QuantReg(Eigen::VectorXd y, Eigen::MatrixXd X, int nObs, int nEqs, Eigen::VectorXd beta, double alpha, Eigen::VectorXd lambda0, int maxIter, double eps);
RcppExport SEXP _bayesELnew_logEL_QuantReg(SEXP ySEXP, SEXP XSEXP, SEXP nObsSEXP, SEXP nEqsSEXP, SEXP betaSEXP, SEXP alphaSEXP, SEXP lambda0SEXP, SEXP maxIterSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type nObs(nObsSEXP);
    Rcpp::traits::input_parameter< int >::type nEqs(nEqsSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(logEL_QuantReg(y, X, nObs, nEqs, beta, alpha, lambda0, maxIter, eps));
    return rcpp_result_gen;
END_RCPP
}
// PostSample_QuantReg
Eigen::MatrixXd PostSample_QuantReg(Eigen::VectorXd y, Eigen::MatrixXd X, int nObs, int nEqs, double alpha, Eigen::VectorXd lambda0, int nsamples, int nburn, Eigen::VectorXd betaInit, Eigen::VectorXd sigs, int maxIter, double eps);
RcppExport SEXP _bayesELnew_PostSample_QuantReg(SEXP ySEXP, SEXP XSEXP, SEXP nObsSEXP, SEXP nEqsSEXP, SEXP alphaSEXP, SEXP lambda0SEXP, SEXP nsamplesSEXP, SEXP nburnSEXP, SEXP betaInitSEXP, SEXP sigsSEXP, SEXP maxIterSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type nObs(nObsSEXP);
    Rcpp::traits::input_parameter< int >::type nEqs(nEqsSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< int >::type nsamples(nsamplesSEXP);
    Rcpp::traits::input_parameter< int >::type nburn(nburnSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type betaInit(betaInitSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type sigs(sigsSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(PostSample_QuantReg(y, X, nObs, nEqs, alpha, lambda0, nsamples, nburn, betaInit, sigs, maxIter, eps));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bayesELnew_lambdaNR", (DL_FUNC) &_bayesELnew_lambdaNR, 9},
    {"_bayesELnew_logEL_MeanReg", (DL_FUNC) &_bayesELnew_logEL_MeanReg, 8},
    {"_bayesELnew_PostSample_MeanReg", (DL_FUNC) &_bayesELnew_PostSample_MeanReg, 11},
    {"_bayesELnew_logEL_QuantReg", (DL_FUNC) &_bayesELnew_logEL_QuantReg, 9},
    {"_bayesELnew_PostSample_QuantReg", (DL_FUNC) &_bayesELnew_PostSample_QuantReg, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_bayesELnew(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
