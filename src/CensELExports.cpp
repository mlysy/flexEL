/**
 * @file CensELExports.cpp
 * 
 * @brief Rcpp wrappers for CensEL class.
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include "inner_elc_new.h"

//[[Rcpp::depends("RcppEigen")]]

/// Construct a CensEL object.
///
/// Instantiates a `CensEL` object on the C++ side, wraps it in an `Rcpp::XPtr`, and returns the corresponding `externalptr` on the R side.
///
/// @param[in] n_obs    Number of observations.
/// @param[in] n_eqs    Number of estimating equations.
/// @return An `externalptr` pointing to the CensEL object.
///
// [[Rcpp::export]]
SEXP CensEL_ctor(int n_obs, int n_eqs) {
  flexEL::CensEL *CEL = new flexEL::CensEL(n_obs, n_eqs);
  Rcpp::XPtr<flexEL::CensEL> pCEL(CEL, true);
  return pCEL;
}

/// Set the maximum number of Newton-Raphson iterations of the CensEL object.
///
/// @param[in] pCEL      `externalptr` pointer to CensEL object. 
/// @param[in] max_iter  Maximum number of Newton-Raphson iterations.
///
// [[Rcpp::export]]
void CensEL_set_max_iter_nr(SEXP pCEL, int max_iter) {
  Rcpp::XPtr<flexEL::CensEL> CEL(pCEL);
  CEL->set_max_iter_nr(max_iter);
  return;
}

/// Set the maximum number of EM iterations of the CensEL object.
///
/// @param[in] pCEL      `externalptr` pointer to CensEL object. 
/// @param[in] max_iter  Maximum number of Newton-Raphson iterations.
///
// [[Rcpp::export]]
void CensEL_set_max_iter_em(SEXP pCEL, int max_iter) {
  Rcpp::XPtr<flexEL::CensEL> CEL(pCEL);
  CEL->set_max_iter_em(max_iter);
  return;
}

/// Set the rel_tol of the CensEL object.
///
/// @param[in] pCEL      `externalptr` pointer to CensEL object. 
/// @param[in] rel_tol   Relative tolerance for the Newton-Raphson algorithm.
///
// [[Rcpp::export]]
void CensEL_set_rel_tol(SEXP pCEL, int rel_tol) {
  Rcpp::XPtr<flexEL::CensEL> CEL(pCEL);
  CEL->set_rel_tol(rel_tol);
  return;
}

/// Set the abs_tol of the CensEL object.
///
/// @param[in] pCEL      `externalptr` pointer to CensEL object. 
/// @param[in] abs_tol   Absolute tolerance for the EM algorithm.
///
// [[Rcpp::export]]
void CensEL_set_abs_tol(SEXP pCEL, int abs_tol) {
  Rcpp::XPtr<flexEL::CensEL> CEL(pCEL);
  CEL->set_rel_tol(abs_tol);
  return;
}

/// Set the supp_adj of the CensEL object.
///
/// @param[in] pCEL        `externalptr` pointer to CensEL object. 
/// @param[in] supp_adj    Whether or not to enable support adjustment.
/// @param[in] a Support   adjustment factor. Defaults to `max(1.0, log(n_obs)/2)`.
///
// [[Rcpp::export]]
void CensEL_set_supp_adj(SEXP pCEL, 
                        bool supp_adj, 
                        Rcpp::Nullable<Rcpp::NumericVector> a_ = R_NilValue) {
  // TODO: is there a better way to handle optional scalar argument?
  Rcpp::XPtr<flexEL::CensEL> CEL(pCEL);
  if (a_.isNull()) {
    CEL->set_supp_adj(supp_adj);
  } else{
    Rcpp::NumericVector a(a_);
    CEL->set_supp_adj(supp_adj, a[0]);
  }
  return;
}



