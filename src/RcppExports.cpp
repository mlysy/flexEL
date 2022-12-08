// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// CensEL_ctor
SEXP CensEL_ctor(int n_obs, int n_eqs);
RcppExport SEXP _flexEL_CensEL_ctor(SEXP n_obsSEXP, SEXP n_eqsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n_obs(n_obsSEXP);
    Rcpp::traits::input_parameter< int >::type n_eqs(n_eqsSEXP);
    rcpp_result_gen = Rcpp::wrap(CensEL_ctor(n_obs, n_eqs));
    return rcpp_result_gen;
END_RCPP
}
// CensEL_get_n_obs
int CensEL_get_n_obs(SEXP pCEL);
RcppExport SEXP _flexEL_CensEL_get_n_obs(SEXP pCELSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pCEL(pCELSEXP);
    rcpp_result_gen = Rcpp::wrap(CensEL_get_n_obs(pCEL));
    return rcpp_result_gen;
END_RCPP
}
// CensEL_get_n_eqs
int CensEL_get_n_eqs(SEXP pCEL);
RcppExport SEXP _flexEL_CensEL_get_n_eqs(SEXP pCELSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pCEL(pCELSEXP);
    rcpp_result_gen = Rcpp::wrap(CensEL_get_n_eqs(pCEL));
    return rcpp_result_gen;
END_RCPP
}
// CensEL_expected_weights
Eigen::VectorXd CensEL_expected_weights(SEXP pCEL, Eigen::VectorXd delta, Eigen::VectorXd epsilon, Eigen::VectorXd omega);
RcppExport SEXP _flexEL_CensEL_expected_weights(SEXP pCELSEXP, SEXP deltaSEXP, SEXP epsilonSEXP, SEXP omegaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pCEL(pCELSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type omega(omegaSEXP);
    rcpp_result_gen = Rcpp::wrap(CensEL_expected_weights(pCEL, delta, epsilon, omega));
    return rcpp_result_gen;
END_RCPP
}
// CensEL_set_smooth
void CensEL_set_smooth(SEXP pCEL, double smooth_s);
RcppExport SEXP _flexEL_CensEL_set_smooth(SEXP pCELSEXP, SEXP smooth_sSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pCEL(pCELSEXP);
    Rcpp::traits::input_parameter< double >::type smooth_s(smooth_sSEXP);
    CensEL_set_smooth(pCEL, smooth_s);
    return R_NilValue;
END_RCPP
}
// CensEL_set_max_iter
void CensEL_set_max_iter(SEXP pCEL, int max_iter);
RcppExport SEXP _flexEL_CensEL_set_max_iter(SEXP pCELSEXP, SEXP max_iterSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pCEL(pCELSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    CensEL_set_max_iter(pCEL, max_iter);
    return R_NilValue;
END_RCPP
}
// CensEL_set_abs_tol
void CensEL_set_abs_tol(SEXP pCEL, double abs_tol);
RcppExport SEXP _flexEL_CensEL_set_abs_tol(SEXP pCELSEXP, SEXP abs_tolSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pCEL(pCELSEXP);
    Rcpp::traits::input_parameter< double >::type abs_tol(abs_tolSEXP);
    CensEL_set_abs_tol(pCEL, abs_tol);
    return R_NilValue;
END_RCPP
}
// GenEL_ctor
SEXP GenEL_ctor(int n_obs, int n_eqs);
RcppExport SEXP _flexEL_GenEL_ctor(SEXP n_obsSEXP, SEXP n_eqsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n_obs(n_obsSEXP);
    Rcpp::traits::input_parameter< int >::type n_eqs(n_eqsSEXP);
    rcpp_result_gen = Rcpp::wrap(GenEL_ctor(n_obs, n_eqs));
    return rcpp_result_gen;
END_RCPP
}
// GenEL_set_max_iter
void GenEL_set_max_iter(SEXP pGEL, int max_iter);
RcppExport SEXP _flexEL_GenEL_set_max_iter(SEXP pGELSEXP, SEXP max_iterSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pGEL(pGELSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    GenEL_set_max_iter(pGEL, max_iter);
    return R_NilValue;
END_RCPP
}
// GenEL_set_rel_tol
void GenEL_set_rel_tol(SEXP pGEL, double rel_tol);
RcppExport SEXP _flexEL_GenEL_set_rel_tol(SEXP pGELSEXP, SEXP rel_tolSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pGEL(pGELSEXP);
    Rcpp::traits::input_parameter< double >::type rel_tol(rel_tolSEXP);
    GenEL_set_rel_tol(pGEL, rel_tol);
    return R_NilValue;
END_RCPP
}
// GenEL_set_supp_adj
void GenEL_set_supp_adj(SEXP pGEL, bool supp_adj, Rcpp::Nullable<Rcpp::NumericVector> a_, Rcpp::Nullable<Rcpp::NumericVector> weight_adj_);
RcppExport SEXP _flexEL_GenEL_set_supp_adj(SEXP pGELSEXP, SEXP supp_adjSEXP, SEXP a_SEXP, SEXP weight_adj_SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pGEL(pGELSEXP);
    Rcpp::traits::input_parameter< bool >::type supp_adj(supp_adjSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type a_(a_SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type weight_adj_(weight_adj_SEXP);
    GenEL_set_supp_adj(pGEL, supp_adj, a_, weight_adj_);
    return R_NilValue;
END_RCPP
}
// GenEL_set_lambda0
void GenEL_set_lambda0(SEXP pGEL, Eigen::VectorXd lambda0);
RcppExport SEXP _flexEL_GenEL_set_lambda0(SEXP pGELSEXP, SEXP lambda0SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pGEL(pGELSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type lambda0(lambda0SEXP);
    GenEL_set_lambda0(pGEL, lambda0);
    return R_NilValue;
END_RCPP
}
// GenEL_get_n_obs
int GenEL_get_n_obs(SEXP pGEL);
RcppExport SEXP _flexEL_GenEL_get_n_obs(SEXP pGELSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pGEL(pGELSEXP);
    rcpp_result_gen = Rcpp::wrap(GenEL_get_n_obs(pGEL));
    return rcpp_result_gen;
END_RCPP
}
// GenEL_get_n_eqs
int GenEL_get_n_eqs(SEXP pGEL);
RcppExport SEXP _flexEL_GenEL_get_n_eqs(SEXP pGELSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pGEL(pGELSEXP);
    rcpp_result_gen = Rcpp::wrap(GenEL_get_n_eqs(pGEL));
    return rcpp_result_gen;
END_RCPP
}
// GenEL_get_supp_adj
bool GenEL_get_supp_adj(SEXP pGEL);
RcppExport SEXP _flexEL_GenEL_get_supp_adj(SEXP pGELSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pGEL(pGELSEXP);
    rcpp_result_gen = Rcpp::wrap(GenEL_get_supp_adj(pGEL));
    return rcpp_result_gen;
END_RCPP
}
// GenEL_get_diag
Rcpp::List GenEL_get_diag(SEXP pGEL);
RcppExport SEXP _flexEL_GenEL_get_diag(SEXP pGELSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pGEL(pGELSEXP);
    rcpp_result_gen = Rcpp::wrap(GenEL_get_diag(pGEL));
    return rcpp_result_gen;
END_RCPP
}
// GenEL_lambda_nr
Eigen::VectorXd GenEL_lambda_nr(SEXP pGEL, Eigen::MatrixXd G, Eigen::VectorXd weights, bool check_conv);
RcppExport SEXP _flexEL_GenEL_lambda_nr(SEXP pGELSEXP, SEXP GSEXP, SEXP weightsSEXP, SEXP check_convSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pGEL(pGELSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type G(GSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< bool >::type check_conv(check_convSEXP);
    rcpp_result_gen = Rcpp::wrap(GenEL_lambda_nr(pGEL, G, weights, check_conv));
    return rcpp_result_gen;
END_RCPP
}
// GenEL_omega_hat
Eigen::VectorXd GenEL_omega_hat(SEXP pGEL, Eigen::VectorXd lambda, Eigen::MatrixXd G, Eigen::VectorXd weights);
RcppExport SEXP _flexEL_GenEL_omega_hat(SEXP pGELSEXP, SEXP lambdaSEXP, SEXP GSEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pGEL(pGELSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type G(GSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(GenEL_omega_hat(pGEL, lambda, G, weights));
    return rcpp_result_gen;
END_RCPP
}
// GenEL_logel
double GenEL_logel(SEXP pGEL, Eigen::MatrixXd G, bool check_conv);
RcppExport SEXP _flexEL_GenEL_logel(SEXP pGELSEXP, SEXP GSEXP, SEXP check_convSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pGEL(pGELSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type G(GSEXP);
    Rcpp::traits::input_parameter< bool >::type check_conv(check_convSEXP);
    rcpp_result_gen = Rcpp::wrap(GenEL_logel(pGEL, G, check_conv));
    return rcpp_result_gen;
END_RCPP
}
// GenEL_weighted_logel
double GenEL_weighted_logel(SEXP pGEL, Eigen::MatrixXd G, Eigen::VectorXd weights, bool check_conv);
RcppExport SEXP _flexEL_GenEL_weighted_logel(SEXP pGELSEXP, SEXP GSEXP, SEXP weightsSEXP, SEXP check_convSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pGEL(pGELSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type G(GSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< bool >::type check_conv(check_convSEXP);
    rcpp_result_gen = Rcpp::wrap(GenEL_weighted_logel(pGEL, G, weights, check_conv));
    return rcpp_result_gen;
END_RCPP
}
// GenEL_logel_grad
Rcpp::List GenEL_logel_grad(SEXP pGEL, Eigen::MatrixXd G, bool check_conv);
RcppExport SEXP _flexEL_GenEL_logel_grad(SEXP pGELSEXP, SEXP GSEXP, SEXP check_convSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pGEL(pGELSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type G(GSEXP);
    Rcpp::traits::input_parameter< bool >::type check_conv(check_convSEXP);
    rcpp_result_gen = Rcpp::wrap(GenEL_logel_grad(pGEL, G, check_conv));
    return rcpp_result_gen;
END_RCPP
}
// GenEL_weighted_logel_grad
Rcpp::List GenEL_weighted_logel_grad(SEXP pGEL, Eigen::MatrixXd G, Eigen::VectorXd weights, bool check_conv);
RcppExport SEXP _flexEL_GenEL_weighted_logel_grad(SEXP pGELSEXP, SEXP GSEXP, SEXP weightsSEXP, SEXP check_convSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pGEL(pGELSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type G(GSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< bool >::type check_conv(check_convSEXP);
    rcpp_result_gen = Rcpp::wrap(GenEL_weighted_logel_grad(pGEL, G, weights, check_conv));
    return rcpp_result_gen;
END_RCPP
}
// adjG
Eigen::MatrixXd adjG(Eigen::MatrixXd G, double a);
RcppExport SEXP _flexEL_adjG(SEXP GSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type G(GSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(adjG(G, a));
    return rcpp_result_gen;
END_RCPP
}
// smooth_indicator
Eigen::VectorXd smooth_indicator(double eps1, Eigen::VectorXd eps2, double s);
RcppExport SEXP _flexEL_smooth_indicator(SEXP eps1SEXP, SEXP eps2SEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type eps1(eps1SEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type eps2(eps2SEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(smooth_indicator(eps1, eps2, s));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_flexEL_CensEL_ctor", (DL_FUNC) &_flexEL_CensEL_ctor, 2},
    {"_flexEL_CensEL_get_n_obs", (DL_FUNC) &_flexEL_CensEL_get_n_obs, 1},
    {"_flexEL_CensEL_get_n_eqs", (DL_FUNC) &_flexEL_CensEL_get_n_eqs, 1},
    {"_flexEL_CensEL_expected_weights", (DL_FUNC) &_flexEL_CensEL_expected_weights, 4},
    {"_flexEL_CensEL_set_smooth", (DL_FUNC) &_flexEL_CensEL_set_smooth, 2},
    {"_flexEL_CensEL_set_max_iter", (DL_FUNC) &_flexEL_CensEL_set_max_iter, 2},
    {"_flexEL_CensEL_set_abs_tol", (DL_FUNC) &_flexEL_CensEL_set_abs_tol, 2},
    {"_flexEL_GenEL_ctor", (DL_FUNC) &_flexEL_GenEL_ctor, 2},
    {"_flexEL_GenEL_set_max_iter", (DL_FUNC) &_flexEL_GenEL_set_max_iter, 2},
    {"_flexEL_GenEL_set_rel_tol", (DL_FUNC) &_flexEL_GenEL_set_rel_tol, 2},
    {"_flexEL_GenEL_set_supp_adj", (DL_FUNC) &_flexEL_GenEL_set_supp_adj, 4},
    {"_flexEL_GenEL_set_lambda0", (DL_FUNC) &_flexEL_GenEL_set_lambda0, 2},
    {"_flexEL_GenEL_get_n_obs", (DL_FUNC) &_flexEL_GenEL_get_n_obs, 1},
    {"_flexEL_GenEL_get_n_eqs", (DL_FUNC) &_flexEL_GenEL_get_n_eqs, 1},
    {"_flexEL_GenEL_get_supp_adj", (DL_FUNC) &_flexEL_GenEL_get_supp_adj, 1},
    {"_flexEL_GenEL_get_diag", (DL_FUNC) &_flexEL_GenEL_get_diag, 1},
    {"_flexEL_GenEL_lambda_nr", (DL_FUNC) &_flexEL_GenEL_lambda_nr, 4},
    {"_flexEL_GenEL_omega_hat", (DL_FUNC) &_flexEL_GenEL_omega_hat, 4},
    {"_flexEL_GenEL_logel", (DL_FUNC) &_flexEL_GenEL_logel, 3},
    {"_flexEL_GenEL_weighted_logel", (DL_FUNC) &_flexEL_GenEL_weighted_logel, 4},
    {"_flexEL_GenEL_logel_grad", (DL_FUNC) &_flexEL_GenEL_logel_grad, 3},
    {"_flexEL_GenEL_weighted_logel_grad", (DL_FUNC) &_flexEL_GenEL_weighted_logel_grad, 4},
    {"_flexEL_adjG", (DL_FUNC) &_flexEL_adjG, 2},
    {"_flexEL_smooth_indicator", (DL_FUNC) &_flexEL_smooth_indicator, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_flexEL(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
