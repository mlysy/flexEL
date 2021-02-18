/**
 * @file inner_el_exports.cpp
 * 
 * @brief Export InnerEL functions to R.
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include "inner_el.h"

//[[Rcpp::depends("RcppEigen")]]

using namespace Eigen;
using namespace Rcpp;

/**
 * @brief Calculate the solution of the dual problem of maximum log EL problem.
 * 
 * @param[in] G          A matrix of dimension `n_eqs x n_obs`.
 * @param[in] max_iter   A positive integer controlling the maximum number of iterations.
 * @param[in] rel_tol    A small positive number controlling accuracy at convergence.
 * @param[in] support    A boolean indicating whether to conduct support correction or not.
 * @param[in] verbose    A boolean indicating whether to print out number of iterations and maximum error 
 *                         at the end of the Newton-Raphson algorithm.
 */
// [[Rcpp::export(".LambdaNR")]]
Eigen::VectorXd LambdaNR(Eigen::MatrixXd G, 
                         Eigen::VectorXd lambda0,
                         int max_iter, 
                         double rel_tol, 
                         bool support, 
                         bool verbose) {
  int n_obs = G.cols();
  int n_eqs = G.rows();
  flexEL::GenEL IL(n_obs,n_eqs);
  IL.set_lambda0(lambda0);
  IL.set_max_iter(max_iter);
  IL.set_rel_tol(rel_tol);
  IL.set_supp_adj(support);
  // IL.set_opts(max_iter,rel_tol,support);
  // IL.set_G(G); // assign the given G
  
  // initialize variables for output here 
  int n_iter;
  double max_err;
  bool not_conv;
  // IL.LambdaNR(n_iter, max_err);
  // VectorXd lambda = IL.get_lambda(); // output (could be not converged)
  VectorXd lambda(n_eqs);
  VectorXd norm_weights = VectorXd::Constant(n_obs+support, 1.0/(n_obs+support));
  // VectorXd norm_weights = VectorXd::Constant(n_obs, 1.0/n_obs);
  IL.lambda_nr(lambda, G, norm_weights);
  IL.get_diag(n_iter, max_err);
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

/**
 * @brief Calculate the probability vector base on the given G matrix.
 * 
 * @param[in] G          A matrix of dimension <code>n_eqs x n_obs</code>.
 * @param[in] max_iter   A positive integer controlling the maximum number of iterations.
 * @param[in] rel_tol    A small positive number controlling accuracy at convergence.
 * @param[in] support    A boolean indicating whether to conduct support correction or not.
 * @param[in] verbose    A boolean indicating whether to print out number of iterations and maximum error.
 *                         at the end of the Newton-Raphson algorithm.
 */
// [[Rcpp::export(".OmegaHat")]]
Eigen::VectorXd OmegaHat(Eigen::MatrixXd G, 
                         Eigen::VectorXd lambda0,
                         int max_iter, 
                         double rel_tol, 
                         bool support, 
                         bool verbose) {
  
  int n_obs = G.cols();
  int n_eqs = G.rows();
  flexEL::GenEL IL(n_obs,n_eqs);
  // IL.set_opts(max_iter,rel_tol,support);
  // IL.set_G(G); // assign the given G
  IL.set_lambda0(lambda0);
  IL.set_max_iter(max_iter);
  IL.set_rel_tol(rel_tol);
  IL.set_supp_adj(support);
  
  // initialize variables for output here 
  int n_iter;
  double max_err;
  bool not_conv;
  // IL.LambdaNR(n_iter, max_err);
  VectorXd lambda(n_eqs);
  VectorXd omega(n_obs + support);
  VectorXd norm_weights = VectorXd::Constant(n_obs+support, 1.0/(n_obs+support));
  IL.lambda_nr(lambda, G, norm_weights);
  IL.get_diag(n_iter, max_err);
  // check convergence
  not_conv = (n_iter == max_iter) && (max_err > rel_tol);
  if(verbose) {
    Rprintf("n_iter = %i, max_err = %f\n", n_iter, max_err);
  }
  
  // int n_obs = G.cols();
  // int n_eqs = G.rows();
  // flexEL::InnerEL IL(n_obs, n_eqs);
  // IL.set_opts(support);
  // IL.set_G(G); // assign the given G
  // IL.set_lambda(lambda); 
  // IL.EvalOmegas(); // calculate omegas
  // VectorXd omegasnew = IL.get_omegas(); // get omegas
  IL.omega_hat(omega, lambda, G, norm_weights);
  return omega;
}

/**
 * @brief Calculate the log empirical likelihood base on the given probability vector.
 * 
 * @param[in] omegas    A numeric probability vector.
 * @param[in] support   A boolean indicating whether to conduct support correction or not.
 */
// [[Rcpp::export(".LogEL")]]
double LogEL(Eigen::VectorXd omegas, bool support) {
  int n_obs = omegas.size() - support;
  flexEL::GenEL IL(n_obs,1);
  // IL.set_opts(support);
  // IL.set_omegas(omegas); 
  // double log_el = IL.LogEL();
  VectorXd norm_weights = VectorXd::Constant(n_obs+support, 1.0/(n_obs+support));
  double sum_weights = double(n_obs + support);
  double log_el = IL.logel_omega(omegas, norm_weights, sum_weights);
  return log_el; 
}

/**
 * @brief Calculate the probability vector, log EL, and the derivative of log EL w.r.t. G evaluated at G.
 * 
 * @param[in] G          A matrix of dimension <code>n_eqs x n_obs</code>.
 * @param[in] max_iter   A positive integer controlling the maximum number of iterations.
 * @param[in] rel_tol    A small positive number controlling accuracy at convergence.
 * @param[in] support    A boolean indicating whether to conduct support correction or not.
 * @param[in] verbose    A boolean indicating whether to print out number of iterations and maximum error.
 *                         at the end of the Newton-Raphson algorithm.
 */
// [[Rcpp::export(".LogELGrad")]]
List LogELGrad(Eigen::MatrixXd G, 
               Eigen::VectorXd lambda0,
               int max_iter, 
               double rel_tol, 
               bool support = false, 
               bool verbose = false) {
  
  int n_obs = G.cols();
  int n_eqs = G.rows();
  // flexEL::InnerEL IL(n_obs,n_eqs);
  // IL.set_opts(max_iter, rel_tol, support);
  // IL.set_G(G); // assign the given G
  flexEL::GenEL IL(n_obs,n_eqs);
  IL.set_lambda0(lambda0);
  IL.set_max_iter(max_iter);
  IL.set_rel_tol(rel_tol);
  IL.set_supp_adj(support);

  int n_iter;
  double max_err;
  // bool not_conv;
  double logel;
  VectorXd lambda(n_eqs);
  VectorXd omega(n_obs + support);
  VectorXd norm_weights = VectorXd::Constant(n_obs+support, 1.0/(n_obs+support));
  double sum_weights = double(n_obs+support);
  MatrixXd dldG(n_eqs, n_obs);
  
  // IL.LambdaNR(n_iter, max_err);
  IL.lambda_nr(lambda, G, norm_weights);
  IL.get_diag(n_iter, max_err);
  if(verbose) {
    Rprintf("n_iter = %i, max_err = %f\n", n_iter, max_err);
  }
  // not_conv = (n_iter == max_iter) && (max_err > rel_tol); // check convergence
  // if (not_conv) {
  //   
  // }
  
  // IL.EvalOmegas();
  // IL.LogELGrad(logel, dldG);
  IL.omega_hat(omega, lambda, G, norm_weights);
  logel = IL.logel_omega(omega, norm_weights, sum_weights);
  IL.logel_grad(dldG, omega, lambda, sum_weights);
  
  return List::create(Named("logel") = logel,
                      Named("dldG") = dldG.transpose(),
                      Named("omega") = omega);
}

