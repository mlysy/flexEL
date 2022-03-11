/// @file example.cpp

// [[Rcpp::depends("flexEL")]]
// [[Rcpp::depends("RcppEigen")]]

#include <Rcpp.h>
#include <RcppEigen.h>
#include <flexEL/gen_el.h>
#include <flexEL/mean_reg_model.h>

/// Empirical likelihood function for the mean regression model.
///
/// @param[in] beta Vector of `n_cov` regression coefficients.
/// @param[in] X Covariate matrix of size `n_obs x n_cov`.
/// @param[in] y Response vector of length `n_obs`.
///
/// @return Scalar value of the log empirical likelihood function.
// [[Rcpp::export]]
double example_logel(Eigen::VectorXd beta,
                     Eigen::MatrixXd X,
                     Eigen::VectorXd y) {
  flexEL::MeanRegModel MR(y, X.transpose());
  int n_obs = MR.get_n_obs();
  int n_eqs = MR.get_n_eqs();
  MatrixXd G = MatrixXd::Zero(MR.get_n_eqs(), MR.get_n_obs());
  MR.eval_G(G, beta);
  flexEL::GenEL GEL(n_obs, n_eqs);
  GEL.set_supp_adj(true); // turn on support correction
  double logel = GEL.logel(G);
  if(!GEL.has_converged_nr()) {
    // convergence check
    logel = -std::numeric_limits<double>::infinity();
  }
  return(logel);
}
