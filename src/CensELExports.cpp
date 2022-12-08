/**
 * @file CensELExports.cpp
 * 
 * @brief Rcpp wrappers for CensEL class.
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include "cens_el.h"

//[[Rcpp::depends("RcppEigen")]]

/// Construct a CensEL object.
///
/// Instantiates a CensEL object on the C++ side, wraps it in an `Rcpp::XPtr`, and returns the corresponding `externalptr` on the R side.
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

/// Getter for `n_obs`.
///
/// @param[in] pCEL `externalptr` pointer to CensEL object.
///
// [[Rcpp::export]]
int CensEL_get_n_obs(SEXP pCEL) {
  Rcpp::XPtr<flexEL::CensEL> CEL(pCEL);
  int n_obs = CEL->get_n_obs();
  return n_obs;
}

/// Getter for `n_eqs`.
///
/// @param[in] pCEL `externalptr` pointer to CensEL object.
/// 
// [[Rcpp::export]]
int CensEL_get_n_eqs(SEXP pCEL) {
  Rcpp::XPtr<flexEL::CensEL> CEL(pCEL);
  int n_eqs = CEL->get_n_eqs();
  return n_eqs;
}


/// Calculate the expected weights in the EM algorithm.
///
/// @param[in] pCEL `externalptr` pointer to CensEL object.
/// @param[in] delta Censoring vector of length `n_obs` or `n_obs+1`.
/// @param[in] epsilon Residual vector of same length as `delta`.
/// @param[in] omega Probability vector of same length as `delta`.
///
/// @return Weight vector of same length as `delta`.
///
// [[Rcpp::export]]
Eigen::VectorXd CensEL_expected_weights(SEXP pCEL,
					Eigen::VectorXd delta,
					Eigen::VectorXd epsilon,
					Eigen::VectorXd omega) {
  Rcpp::XPtr<flexEL::CensEL> CEL(pCEL);
  int n_obs2 = delta.size(); // input checking is done on the R side
  Eigen::VectorXd weights(n_obs2);
  // Eigen::Ref<const Eigen::VectorXd> delta_eff = CEL->supp_delta(delta);
  // Eigen::Ref<const Eigen::VectorXd> epsilon_eff = CEL->supp_epsilon(epsilon);
  CEL->expected_weights(weights, delta, epsilon, omega);
  return weights;
}

/// Set the smoothing parameter of the CensEL object.
///
/// @param[in] pCEL     `externalptr` pointer to CensEL object. 
/// @param[in] smooth_s Smoothing parameter (a positive scalar). Defaults to 10.
///
// [[Rcpp::export]]
void CensEL_set_smooth(SEXP pCEL, double smooth_s) {
  Rcpp::XPtr<flexEL::CensEL> CEL(pCEL);
  CEL->set_smooth(smooth_s);
  return;
}

/// Set the maximum number of EM iterations.
///
/// @param[in] pCEL     `externalptr` pointer to CensEL object. 
/// @param[in] max_iter Maximum number of EM iterations.
///
// [[Rcpp::export]]
void CensEL_set_max_iter(SEXP pCEL, int max_iter) {
  Rcpp::XPtr<flexEL::CensEL> CEL(pCEL);
  CEL->set_max_iter(max_iter);
  return;
}

/// Set the absolute tolerance for terminating the EM algorithm.
///
/// @param[in] pCEL     `externalptr` pointer to CensEL object. 
/// @param[in] abs_tol Absolute tolerance.
///
// [[Rcpp::export]]
void CensEL_set_abs_tol(SEXP pCEL, double abs_tol) {
  Rcpp::XPtr<flexEL::CensEL> CEL(pCEL);
  CEL->set_abs_tol(abs_tol);
  return;
}

