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

/// Getter for n_obs.
///
/// @param[in] pGEL   `externalptr` pointer to CensEL object. 
// [[Rcpp::export]]
int CensEL_get_n_obs(SEXP pGEL) {
  Rcpp::XPtr<flexEL::CensEL> GEL(pGEL);
  int n_obs = GEL->get_n_obs();
  return n_obs;
}

/// Getter for n_eqs.
///
/// @param[in] pGEL   `externalptr` pointer to CensEL object. 
// [[Rcpp::export]]
int CensEL_get_n_eqs(SEXP pGEL) {
  Rcpp::XPtr<flexEL::CensEL> GEL(pGEL);
  int n_eqs = GEL->get_n_eqs();
  return n_eqs;
}

/// Getter for supp_adj.
///
/// @param[in] pGEL   `externalptr` pointer to CensEL object. 
// [[Rcpp::export]]
bool CensEL_get_supp_adj(SEXP pGEL) {
  Rcpp::XPtr<flexEL::CensEL> GEL(pGEL);
  bool supp_adj = GEL->get_supp_adj();
  return supp_adj;
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
  Rcpp::XPtr<flexEL::CensEL> CEL(pCEL);
  // TODO: is there a better way to handle optional scalar argument?
  if (a_.isNull()) {
    CEL->set_supp_adj(supp_adj);
  } else{
    Rcpp::NumericVector a(a_);
    CEL->set_supp_adj(supp_adj, a[0]);
  }
  return;
}

/// ...
///
/// @param[in] pCEL     `externalptr` pointer to CensEL object. 
/// @param[in] G        Moment matrix of size `n_eqs x n_obs` or `n_eqs x (n_obs + supp_adj)`.  If `supp_adj = false`, the former is required.  If `supp_adj = true` and the former is provided, support adjustment is performed.  If `supp_adj = true` and `G.cols() == n_obs + 1`, assumes that support has already been corrected. 
/// @param[in] delta    Vector of censoring indicator of length `n_obs`.
/// @param[in] epsilon  Vector of residuals of length `n_obs`.
///
// [[Rcpp::export]]
Eigen::VectorXd CensEL_eval_weights(SEXP pCEL,
                                    Eigen::VectorXd delta,
                                    Eigen::VectorXd epsilon,
                                    Eigen::VectorXd omega) {
  Rcpp::XPtr<flexEL::CensEL> CEL(pCEL);
  int n_obs = CEL->get_n_obs(); // G.cols() should be the same value, check in R side
  // std::cout << "CensEL_omega_hat: n_obs = " << n_obs << std::endl;
  bool supp_adj = CEL->get_supp_adj();
  // std::cout << "CensEL_omega_hat: supp_adj = " << supp_adj << std::endl;
  Eigen::VectorXd weights(n_obs + supp_adj);
  CEL->eval_weights(weights, delta, epsilon, omega);
  // CEL->omega_hat(omega, G, delta, epsilon);
  return weights;
}

/// ...
///
/// @param[in] pCEL     `externalptr` pointer to CensEL object.
/// @param[in] G        Moment matrix of size `n_eqs x n_obs` or `n_eqs x (n_obs + supp_adj)`.  If `supp_adj = false`, the former is required.  If `supp_adj = true` and the former is provided, support adjustment is performed.  If `supp_adj = true` and `G.cols() == n_obs + 1`, assumes that support has already been corrected.
/// @param[in] delta    Vector of censoring indicator of length `n_obs`.
/// @param[in] epsilon  Vector of residuals of length `n_obs`.
///
// [[Rcpp::export]]
Eigen::VectorXd CensEL_omega_hat(SEXP pCEL,
                      Eigen::MatrixXd G,
                      Eigen::VectorXd delta,
                      Eigen::VectorXd epsilon) {
  Rcpp::XPtr<flexEL::CensEL> CEL(pCEL);
  int n_obs = CEL->get_n_obs(); // G.cols() should be the same value, check in R side
  // std::cout << "CensEL_omega_hat: n_obs = " << n_obs << std::endl;
  bool supp_adj = CEL->get_supp_adj();
  // std::cout << "CensEL_omega_hat: supp_adj = " << supp_adj << std::endl;
  Eigen::VectorXd omega =  Eigen::VectorXd::Constant(n_obs+supp_adj, 1.0/(n_obs+supp_adj));
  CEL->omega_hat(omega, G, delta, epsilon);
  return omega;
}

