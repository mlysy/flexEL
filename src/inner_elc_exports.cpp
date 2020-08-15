/**
 * @file inner_elc_exports.cpp
 * 
 * @brief Export InnerELC functions to R.
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include "inner_elc.h"

//[[Rcpp::depends("RcppEigen")]]

using namespace Rcpp;
using namespace Eigen;

/**
 * @brief Calculate the weights for the weighted log EL.
 * 
 * @param[in] omegas   A numeric probability vector of the same length as \c deltas.
 * @param[in] deltas   A numeric vector of 0 and 1s where 1 indicates fully observed value and 0 indicates right-censored value.
 * @param[in] epsilons A numeric vector of residuals of the same length as \c deltas.
 * @param[in] support  A boolean indicating whether to conduct support correction or not.
 */
// [[Rcpp::export(".EvalWeights")]]
Eigen::VectorXd EvalWeights(Eigen::VectorXd omegas, 
                            Eigen::VectorXd deltas,
                            Eigen::VectorXd epsilons, 
                            bool support) {
  int n_obs = epsilons.size(); 
  flexEL::InnerELC ILC(n_obs,1);
  ILC.set_opts(support);
  ILC.set_deltas(deltas);
  ILC.set_omegas(omegas);
  ILC.set_epsilons(epsilons); 
  ILC.EvalWeights(); 
  Eigen::VectorXd weights = ILC.get_weights(); 
  return(weights);
}

/**
 * @brief Calculate the solution of the dual problem of the weighted maximum log EL problem.
 * 
 * @param[in] G          A matrix of dimension <code>n_eqs x n_obs</code>.
 * @param[in] weights    A numeric vector serving as the weights of the weighted maximum log EL problem.
 * @param[in] max_iter   A positive integer controlling the maximum number of iterations.
 * @param[in] rel_tol    A small positive number controlling accuracy at convergence.
 * @param[in] support    A boolean indicating whether to conduct support correction or not.
 * @param[in] verbose    A boolean indicating whether to print out number of iterations and maximum error.
 */
// [[Rcpp::export(".LambdaNRCens")]]
Eigen::VectorXd LambdaNRCens(Eigen::MatrixXd G, 
                             Eigen::VectorXd weights, 
                             int max_iter, 
                             double rel_tol, 
                             bool support, 
                             bool verbose) {
  int n_obs = G.cols();
  int n_eqs = G.rows();
  flexEL::InnerELC ILC(n_obs,n_eqs);
  ILC.set_opts(max_iter, rel_tol, support);
  ILC.set_G(G); // assign a given G
  ILC.set_weights(weights);

  // initialize variables for output here 
  int n_iter;
  double max_err;
  ILC.LambdaNR(n_iter, max_err);
  VectorXd lambda = ILC.get_lambda(); // output

  // check convergence
  bool not_conv = (n_iter == max_iter) && (max_err > rel_tol);
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
 * @brief Calculate the probability vector given right-censored observation using an EM algorithm.
 * 
 * @param[in] G           A matrix of dimension <code>n_eqs x n_obs</code>.
 * @param[in] omegas_init A numeric probability vector serving as initial value.
 * @param[in] deltas   A numeric vector of 0 and 1s where 1 indicates fully observed value and 0 indicates right-censored value.
 * @param[in] epsilons A numeric vector of residuals of the same length as \c deltas.
 * @param[in] max_iter    A positive integer controlling the maximum number of iterations.
 * @param[in] rel_tol     A small positive number controlling accuracy at convergence of the Newton-Raphson algorithm (tolerance of relative error).
 * @param[in] abs_tol     A small positive number controlling accuracy at convergence of the EM algorithm (tolerance of absolute error).
 * @param[in] support     A boolean indicating whether to conduct support correction or not.
 * @param[in] verbose     A boolean indicating whether to print out number of iterations and maximum error.
 *                          at the end of the Newton-Raphson algorithm.
 */
// [[Rcpp::export(".OmegaHatEM")]]
Eigen::VectorXd OmegaHatEM(Eigen::MatrixXd G, 
                           Eigen::VectorXd omegas_init, 
                           Eigen::VectorXd deltas,
                           Eigen::VectorXd epsilons, 
                           int max_iter, double rel_tol, double abs_tol, 
                           bool support, 
                           bool verbose) {
  int n_obs = G.cols();
  int n_eqs = G.rows();
  flexEL::InnerELC ILC(n_obs,n_eqs); 
  ILC.set_opts(max_iter, rel_tol, abs_tol, support);
  ILC.set_deltas(deltas);
  ILC.set_G(G); // assign a given G
  ILC.set_epsilons(epsilons);
  ILC.set_omegas(omegas_init); // set initial omegas from uncensored omega.hat
  ILC.EvalOmegas();
  VectorXd omegasnew = ILC.get_omegas(); // output
  return omegasnew;
}

/**
 * @brief Calculate the censored log EL.
 * 
 * @param[in] omegas   A numeric probability vector of the same length as \c deltas.
 * @param[in] deltas   A numeric vector of 0 and 1s where 1 indicates fully observed value and 0 indicates right-censored value.
 * @param[in] epsilons A numeric vector of residuals of the same length as \c deltas.
 * @param[in] support  A boolean indicating whether to conduct support correction or not.
 */
// [[Rcpp::export(".LogELCens")]]
double LogELCens(Eigen::VectorXd omegas, 
                 Eigen::VectorXd deltas, 
                 Eigen::VectorXd epsilons, 
                 bool support) {
  int n_obs = omegas.size() - support;
  flexEL::InnerELC ILC(n_obs,1);
  ILC.set_opts(support);
  ILC.set_deltas(deltas);
  ILC.set_epsilons(epsilons); 
  ILC.set_omegas(omegas);
  double logel = ILC.LogEL();
  return logel; 
}

/**
 * @brief Calculate the censored log EL with continuity correction.
 * 
 * @param[in] omegas   A numeric probability vector of the same length as \c deltas.
 * @param[in] deltas   A numeric vector of 0 and 1s where 1 indicates fully observed value and 0 indicates right-censored value.
 * @param[in] epsilons A numeric vector of residuals of the same length as \c deltas.
 * @param[in] sp       A numeric scalar as the tuning parameter for continuity correction.
 * @param[in] support  A boolean indicating whether to conduct support correction or not.
 */
// [[Rcpp::export(".LogELSmooth")]]
double LogELSmooth(Eigen::VectorXd omegas, 
                   Eigen::VectorXd deltas,
                   Eigen::VectorXd epsilons, 
                   double sp, 
                   bool support) {
  int n_obs = omegas.size() - support;
  flexEL::InnerELC ILC(n_obs,1);
  ILC.set_opts(support);
  ILC.set_omegas(omegas);
  ILC.set_epsilons(epsilons);
  ILC.set_deltas(deltas);
  return ILC.LogELSmooth(sp);
}

/**
 * @brief Calculate the weights for the weighted log EL with continuity correction.
 * 
 * @param[in] omegas   A numeric probability vector of the same length as \c deltas.
 * @param[in] deltas   A numeric vector of 0 and 1s where 1 indicates fully observed value and 0 indicates right-censored value.
 * @param[in] epsilons A numeric vector of residuals of the same length as \c deltas.
 * @param[in] sp       A numeric scalar as the tuning parameter for continuity correction.
 * @param[in] support  A boolean indicating whether to conduct support correction or not.
 */
// [[Rcpp::export(".EvalWeightsSmooth")]]
Eigen::VectorXd EvalWeightsSmooth(Eigen::VectorXd omegas,
                                  Eigen::VectorXd deltas, 
                                  Eigen::VectorXd epsilons, 
                                  double sp, 
                                  bool support) {
  int n_obs = epsilons.size();
  flexEL::InnerELC ILC(n_obs,1);
  ILC.set_opts(support);
  ILC.set_omegas(omegas);
  ILC.set_epsilons(epsilons);
  ILC.set_deltas(deltas);
  ILC.EvalWeightsSmooth(sp);
  Eigen::VectorXd weights = ILC.get_weights(); 
  return(weights);
}


/**
 * @brief Calculate the probability vector given right-censored observation using an EM algorithm.
 * 
 * @param[in] G           A matrix of dimension <code>n_eqs x n_obs</code>.
 * @param[in] omegas_init A numeric probability vector serving as initial value.
 * @param[in] deltas      A numeric vector of 0 and 1s where 1 indicates fully observed value and 0 indicates right-censored value.
 * @param[in] epsilons    A numeric vector of residuals of the same length as \c deltas.
 * @param[in] sp          A numeric scalar as the tuning parameter for continuity correction.
 * @param[in] max_iter    A positive integer controlling the maximum number of iterations.
 * @param[in] rel_tol     A small positive number controlling accuracy at convergence of the Newton-Raphson algorithm (tolerance of relative error).
 * @param[in] abs_tol     A small positive number controlling accuracy at convergence of the EM algorithm (tolerance of absolute error).
 * @param[in] support     A boolean indicating whether to conduct support correction or not.
 * @param[in] verbose     A boolean indicating whether to print out number of iterations and maximum error.
 *                          at the end of the Newton-Raphson algorithm.
 */
// [[Rcpp::export(".OmegaHatEMSmooth")]]
Eigen::VectorXd OmegaHatEMSmooth(Eigen::MatrixXd G, 
                                 Eigen::VectorXd omegas_init,
                                 Eigen::VectorXd deltas,
                                 Eigen::VectorXd epsilons, 
                                 double sp, 
                                 int max_iter, double rel_tol, double abs_tol, 
                                 bool support, 
                                 bool verbose) {
  int n_obs = G.cols();
  int n_eqs = G.rows();
  flexEL::InnerELC ILC(n_obs,n_eqs); 
  ILC.set_opts(max_iter, rel_tol, abs_tol, support);
  ILC.set_deltas(deltas);
  ILC.set_G(G); // assign a given G
  ILC.set_epsilons(epsilons); 
  // ILC.setTol(max_iter, rel_tol, abs_tol);
  ILC.set_omegas(omegas_init); // set initial omegas from uncensored omega.hat
  ILC.EvalOmegasSmooth(sp);
  VectorXd omegasnew = ILC.get_omegas(); // output
  return omegasnew; 
}
