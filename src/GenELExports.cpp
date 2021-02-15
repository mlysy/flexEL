/**
 * @file GenELExports.cpp
 * 
 * @brief Rcpp wrappers for GenEL class.
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include "inner_el.h"

//[[Rcpp::depends("RcppEigen")]]

/// Construct a GenEL object.
///
/// Instantiates a `GenEL` object on the C++ side, wraps it in an `Rcpp::XPtr`, and returns the corresponding `externalptr` on the R side.
///
/// @param[in] n_obs    Number of observations.
/// @param[in] n_eqs    Number of estimating equations.
/// @return An `externalptr` pointing to the GenEL object.
///
// [[Rcpp::export]]
SEXP GenEL_ctor(int n_obs, int n_eqs) {
  flexEL::GenEL *GEL = new flexEL::GenEL(n_obs, n_eqs);
  Rcpp::XPtr<flexEL::GenEL> pGEL(GEL, true);
  return pGEL;
}

/// Set the max_iter of the GenEL object.
///
/// @param[in] pGEL `externalptr` pointer to GenEL object. 
/// @param[in] max_iter Maximum number of Newton-Raphson iterations.
///
// [[Rcpp::export]]
void GenEL_set_max_iter(SEXP pGEL, int max_iter) {
  Rcpp::XPtr<flexEL::GenEL> GEL(pGEL);
  GEL->set_max_iter(max_iter);
  return;
}

/// Set the rel_tol of the GenEL object.
///
/// @param[in] pGEL `externalptr` pointer to GenEL object. 
/// @param[in] rel_tol Relative tolerance for the Newton-Raphson algorithm.
///
// [[Rcpp::export]]
void GenEL_set_rel_tol(SEXP pGEL, int rel_tol) {
  Rcpp::XPtr<flexEL::GenEL> GEL(pGEL);
  GEL->set_rel_tol(rel_tol);
  return;
}

/// Set the supp_adj of the GenEL object.
///
/// @param[in] pGEL `externalptr` pointer to GenEL object. 
/// @param[in] supp_adj Whether or not to enable support adjustment.
/// @param[in] a Support adjustment factor. Defaults to `max(1.0, log(n_obs)/2)`.
///
// [[Rcpp::export]]
void GenEL_set_supp_adj(SEXP pGEL, 
                        bool supp_adj, 
                        Rcpp::Nullable<Rcpp::NumericVector> a_ = R_NilValue) {
  // TODO: is there a better way to handle optional scalar argument?
  Rcpp::XPtr<flexEL::GenEL> GEL(pGEL);
  if (a_.isNull()) {
    GEL->set_supp_adj(supp_adj);
  } else{
    Rcpp::NumericVector a(a_);
    GEL->set_supp_adj(supp_adj, a[0]);
  }
  return;
}

/// Set the lambda0 of the GenEL object.
///
/// @param[in] pGEL `externalptr` pointer to GenEL object. 
/// @param[in] lambda0 Initialization vector of size `n_eqs`.
///
// [[Rcpp::export]]
void GenEL_set_lambda0(SEXP pGEL, Eigen::VectorXd lambda0) {
  Rcpp::XPtr<flexEL::GenEL> GEL(pGEL);
  GEL->set_lambda0(lambda0);
  return;
}

/// Solve the dual problem via Newton-Raphson algorithm.
///
/// @param[in] pGEL `externalptr` pointer to GenEL object. 
/// @param[in] G Moment matrix of size `n_eqs x n_obs` or `n_eqs x (n_obs + supp_adj)`.  If `supp_adj = false`, the former is required.  If `supp_adj = true` and the former is provided, support adjustment is performed.  If `supp_adj = true` and `G.cols() == n_obs + 1`, assumes that support has already been corrected. 
///
// [[Rcpp::export]]
Eigen::VectorXd GenEL_lambda_nr(SEXP pGEL, Eigen::MatrixXd G, bool verbose) {
  Rcpp::XPtr<flexEL::GenEL> GEL(pGEL);
  bool supp_adj = GEL->get_supp_adj();
  int n_obs = G.cols();
  int n_eqs = G.rows();
  Eigen::VectorXd lambda(n_eqs);
  Eigen::VectorXd norm_weights = Eigen::VectorXd::Constant(n_obs+supp_adj, 1.0/(n_obs+supp_adj));
  GEL->lambda_nr(lambda, G, norm_weights);
  int n_iter;
  double max_err;
  bool not_conv;
  GEL->get_diag(n_iter, max_err);
  not_conv = (n_iter == GEL->get_max_iter()) && (max_err > GEL->get_rel_tol());
  if(verbose) {
    Rprintf("n_iter = %i, max_err = %f\n", n_iter, max_err);
  }
  if (not_conv) {
    for (int ii=0; ii < lambda.size(); ii++) {
      lambda(ii) = std::numeric_limits<double>::quiet_NaN();
    }
  }
  return lambda;
}

/// Getter for n_obs.
///
/// @param[in] pGEL `externalptr` pointer to GenEL object. 
// [[Rcpp::export]]
int GenEL_get_n_obs(SEXP pGEL) {
  Rcpp::XPtr<flexEL::GenEL> GEL(pGEL);
  int n_obs = GEL->get_n_obs();
  return n_obs;
}

/// Getter for n_eqs.
///
/// @param[in] pGEL `externalptr` pointer to GenEL object. 
// [[Rcpp::export]]
int GenEL_get_n_eqs(SEXP pGEL) {
  Rcpp::XPtr<flexEL::GenEL> GEL(pGEL);
  int n_eqs = GEL->get_n_eqs();
  return n_eqs;
}

/// Calculate the probability vector base on the given G matrix.
/// 
/// @param[in] pGEL `externalptr` pointer to GenEL object. 
/// @param[in] lambda Dual problem vector of size `n_eqs`.  
/// @param[in] G Moment matrix of size `n_eqs x n_obs` or `n_eqs x (n_obs + supp_adj)`.  If `supp_adj = false`, the former is required.  If `supp_adj = true` and the former is provided, support adjustment is performed.  If `supp_adj = true` and `G.cols() == n_obs + 1`, assumes that support has already been corrected. 
// [[Rcpp::export(".OmegaHat")]]
Eigen::VectorXd GenEL_omega_hat(SEXP pGEL, 
                                Eigen::VectorXd lambda,
                                Eigen::MatrixXd G) {
  Rcpp::XPtr<flexEL::GenEL> GEL(pGEL);
  int n_obs = GEL->get_n_obs(); // G.cols() should be the same value, check in R side
  bool supp_adj = GEL->get_supp_adj();
  Eigen::VectorXd omega(n_obs + supp_adj);
  Eigen::VectorXd norm_weights = Eigen::VectorXd::Constant(n_obs+supp_adj, 1.0/(n_obs+supp_adj));
  GEL->omega_hat(omega, lambda, G, norm_weights);
  return omega;
}

/// Calculate the log empirical likelihood base on the given probability vector.
/// 
/// @param[in] pGEL `externalptr` pointer to GenEL object. 
/// @param[in] omega Probability vector of length `n_obs + supp_adj`.
// [[Rcpp::export(".OmegaHat")]]
double GenEL_logel_omega(SEXP pGEL,
                         Eigen::VectorXd omega) {
  Rcpp::XPtr<flexEL::GenEL> GEL(pGEL);
  int n_obs = GEL->get_n_obs(); // length of omega should be the same, check in R side
  bool supp_adj = GEL->get_supp_adj();
  Eigen::VectorXd norm_weights = Eigen::VectorXd::Constant(n_obs+supp_adj, 1.0/(n_obs+supp_adj));
  double sum_weights = double(n_obs + supp_adj);
  double log_el = GEL->logel_omega(omega, norm_weights, sum_weights);
  return log_el;
}


  