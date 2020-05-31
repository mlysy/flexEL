/**
 * @file InnerElExports.cpp
 * 
 * @brief Export InnerEL functions to R.
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include "inner_el.h"

//[[Rcpp::depends("RcppEigen")]]

using namespace Eigen;
using namespace Rcpp;

// [[Rcpp::export(".LambdaNR")]]
Eigen::VectorXd LambdaNR(Eigen::MatrixXd G, int max_iter, double rel_tol, bool support, bool verbose) {
  int n_obs = G.cols();
  int n_eqs = G.rows();
  flexEL::InnerEL IL(n_obs,n_eqs);
  IL.set_opts(max_iter,rel_tol,support);
  IL.set_G(G); // assign the given G
  
  // initialize variables for output here 
  int n_iter;
  double max_err;
  bool not_conv;
  // IL.setTol(max_iter, rel_tol);
  IL.LambdaNR(n_iter, max_err);
  VectorXd lambda = IL.get_lambda(); // output (could be not converged)
  // check convergence
  not_conv = (n_iter == max_iter) && (max_err > rel_tol);
  if(verbose) {
    Rprintf("n_iter = %i, max_err = %f\n", n_iter, max_err);
  }
  
  // fill in NaN if not converged
  if (not_conv) {
    for (int ii=0; ii<lambda.size(); ii++) {
      lambda(ii) = std::numeric_limits<double>::quiet_NaN();
    }
  }
  return lambda;
}

// Eigen::VectorXd y, Eigen::MatrixXd X not needed here since there is no ordering 
// as in the censored case 
// [[Rcpp::export(".OmegaHat")]]
Eigen::VectorXd OmegaHat(Eigen::MatrixXd G, Eigen::VectorXd lambda, bool support) {
  int n_obs = G.cols();
  int n_eqs = G.rows();
  flexEL::InnerEL IL(n_obs,n_eqs);
  IL.set_opts(support);
  IL.set_G(G); // assign the given G
  IL.set_lambda(lambda); 
  IL.EvalOmegas(); // calculate omegas
  VectorXd omegasnew = IL.get_omegas(); // get omegas
  return omegasnew;
}

// Calculate logEL given omegas
// [[Rcpp::export(".LogEL")]]
double LogEL(Eigen::VectorXd omegas, bool support) {
  int n_obs = omegas.size();
  flexEL::InnerEL IL(n_obs,1);
  IL.set_opts(support);
  IL.set_omegas(omegas); 
  double log_el = IL.LogEL();
  return log_el; 
}

