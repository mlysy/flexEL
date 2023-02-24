/**
 * @file utilsExports.cpp
 * 
 * @brief Export functions from `utils.h` for testing purposes.
 */

#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends("RcppEigen")]]
#include <RcppEigen.h>
using namespace Eigen;
#include "utils.h"

/// Calculate the adjusted `G` matrix.
///
/// @param[in] G A numeric matrix of size `n_eqs x n_obs`.
/// @param[in] a Scalar tuning parameter for the adjustment.
/// @return The adjusted `G` matrix of size `n_eqs x (n_obs + 1)`.
///
// [[Rcpp::export(".adjG")]]
Eigen::MatrixXd adjG(Eigen::MatrixXd G, double a) {
  int n_obs = G.cols();
  int n_eqs = G.rows();
  Eigen::MatrixXd aG = Eigen::MatrixXd::Zero(n_eqs, n_obs+1);
  aG.leftCols(n_obs) = G;
  flexEL::adj_G<double>(aG, a);
  return aG;
}

/// Smoothed indicator function.
///
/// Calculates `smooth_ind = 1/(1 + exp(s * (eps1 - eps2)))`, where `eps1` is a scalar and `eps2` is a vector.
///
/// @param eps1[in] First residual.
/// @param eps2[in] Second residual.
/// @param s[in] Smoothing parameter.
/// @return Smoothed indicator vector the same size as `eps2`.
///
// [[Rcpp::export]]
Eigen::VectorXd smooth_indicator(double eps1,
				 Eigen::VectorXd eps2,
				 double s) {
  Eigen::VectorXd smooth_ind = Eigen::VectorXd::Zero(eps2.size());
  flexEL::smooth_indicator<double>(smooth_ind, eps1, eps2, s);
  return smooth_ind;
}

