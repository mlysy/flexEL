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
  flexEL::CensEL<double> *cel = new flexEL::CensEL<double>(n_obs, n_eqs);
  Rcpp::XPtr<flexEL::CensEL<double> > p_cel(cel, true);
  return p_cel;
}

/// Getter for `n_obs`.
///
/// @param[in] p_cel `externalptr` pointer to CensEL object.
///
// [[Rcpp::export]]
int CensEL_get_n_obs(SEXP p_cel) {
  Rcpp::XPtr<flexEL::CensEL<double> > cel(p_cel);
  int n_obs = cel->get_n_obs();
  return n_obs;
}

/// Getter for `n_eqs`.
///
/// @param[in] p_cel `externalptr` pointer to CensEL object.
/// 
// [[Rcpp::export]]
int CensEL_get_n_eqs(SEXP p_cel) {
  Rcpp::XPtr<flexEL::CensEL<double> > cel(p_cel);
  int n_eqs = cel->get_n_eqs();
  return n_eqs;
}


/// Calculate the expected weights in the EM algorithm.
///
/// @param[in] p_cel `externalptr` pointer to CensEL object.
/// @param[in] delta Censoring vector of length `n_obs` or `n_obs+1`.
/// @param[in] epsilon Residual vector of same length as `delta`.
/// @param[in] omega Probability vector of same length as `delta`.
///
/// @return Weight vector of same length as `delta`.
///
// [[Rcpp::export]]
Eigen::VectorXd CensEL_expected_weights(SEXP p_cel,
					Eigen::VectorXd delta,
					Eigen::VectorXd epsilon,
					Eigen::VectorXd omega) {
  Rcpp::XPtr<flexEL::CensEL<double> > cel(p_cel);
  int n_obs2 = delta.size(); // input checking is done on the R side
  Eigen::VectorXd weights(n_obs2);
  // Eigen::Ref<const Eigen::VectorXd> delta_eff = cel->supp_delta(delta);
  // Eigen::Ref<const Eigen::VectorXd> epsilon_eff = cel->supp_epsilon(epsilon);
  cel->expected_weights(weights, delta, epsilon, omega);
  return weights;
}

/// Calculate the profile probability weights.
///
/// @param[in] p_cel `externalptr` pointer to CensEL object.
/// @param[in] G Moment matrix of size `n_eqs x n_obs`.
/// @param[in] delta Censoring vector of length `n_obs`.
/// @param[in] epsilon Residual vector of length `n_obs`.
/// @param[in] check_conv If `true`, checks whether each NR calculation converged and also that `abs_tol` was reached within the given number of EM iterations.  If not, sets the probability weights to a vector of `NaN`s.
/// @param[in] p_gel `externalptr` pointer to GenEL object.
///
/// @return Probability vector of length `n_obs + supp_adj`.
///
// [[Rcpp::export]]
Eigen::VectorXd CensEL_omega_hat(SEXP p_cel,
				 Eigen::MatrixXd G,
				 Eigen::VectorXd delta,
				 Eigen::VectorXd epsilon,
				 bool check_conv,
				 SEXP p_gel) {
  Rcpp::XPtr<flexEL::CensEL<double> > cel(p_cel);
  Rcpp::XPtr<flexEL::GenEL<double> > gel(p_gel);
  int n_obs = gel->get_n_obs();
  bool supp_adj = gel->get_supp_adj();
  Eigen::VectorXd omega(n_obs + supp_adj);
  cel->omega_hat(omega, G, delta, epsilon, *gel);
  bool has_conv = check_conv ? cel->has_converged_em() : true;
  if(!has_conv) {
    omega.setConstant(std::numeric_limits<double>::quiet_NaN());
  }
  return omega;
}


/// Set the smoothing parameter of the CensEL object.
///
/// @param[in] p_cel     `externalptr` pointer to CensEL object. 
/// @param[in] smooth_s Smoothing parameter (a positive scalar). Defaults to 10.
///
// [[Rcpp::export]]
void CensEL_set_smooth(SEXP p_cel, double smooth_s) {
  Rcpp::XPtr<flexEL::CensEL<double> > cel(p_cel);
  cel->set_smooth(smooth_s);
  return;
}

/// Set the maximum number of EM iterations.
///
/// @param[in] p_cel     `externalptr` pointer to CensEL object. 
/// @param[in] max_iter Maximum number of EM iterations.
///
// [[Rcpp::export]]
void CensEL_set_max_iter(SEXP p_cel, int max_iter) {
  Rcpp::XPtr<flexEL::CensEL<double> > cel(p_cel);
  cel->set_max_iter(max_iter);
  return;
}

/// Set the absolute tolerance for terminating the EM algorithm.
///
/// @param[in] p_cel     `externalptr` pointer to CensEL object. 
/// @param[in] abs_tol Absolute tolerance.
///
// [[Rcpp::export]]
void CensEL_set_abs_tol(SEXP p_cel, double abs_tol) {
  Rcpp::XPtr<flexEL::CensEL<double> > cel(p_cel);
  cel->set_abs_tol(abs_tol);
  return;
}

