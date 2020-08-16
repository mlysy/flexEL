/**
 * @file mean_reg_model_exports.cpp
 * 
 * @brief Export mean_reg_model functions to R.
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include "mean_reg_model.h"

//[[Rcpp::depends("RcppEigen")]]

// using namespace Rcpp;
using namespace Eigen;


/**
 * @brief Calculate the G matrix for regular mean regression model.
 * 
 * @param[in] y      Responses of length <code>n_obs</code>.
 * @param[in] X      Covariate matrix of dimension <code>nBet</code> x <code>n_obs</code>.
 * @param[in] beta   Coefficient vector of length <code>n_bet_</code> in linear location function.
 */
// [[Rcpp::export(".MeanRegEvalG")]]
Eigen::MatrixXd MeanReg_evalG(Eigen::VectorXd y, 
                              Eigen::MatrixXd X, 
                              Eigen::VectorXd beta) {
  flexEL::MeanRegModel MR(y, X);
  MatrixXd G = MatrixXd::Zero(MR.get_n_eqs(), MR.get_n_obs()); // G in C++ is n_eqs * n_obs
  MR.EvalG(G, beta);
  return(G);
}

/**
 * @brief Calculate the G matrix for location-scale mean regression model.
 * 
 * @param[in] y      Responses of length <code>n_obs</code>.
 * @param[in] X      Covariate matrix of dimension <code>nBet</code> x <code>n_obs</code>.
 * @param[in] Z      Covariate matrix of dimension <code>nGam</code> x <code>n_obs</code>.
 * @param[in] beta   Coefficient vector of length <code>nBet</code> in linear location function.
 * @param[in] gamma  Coefficient vector of length <code>nGam</code> in exponential scale function.
 * @param[in] sig2   Scale parameter in scale function.
 */
// [[Rcpp::export(".MeanRegLSEvalG")]]
Eigen::MatrixXd MeanRegLS_EvalG(Eigen::VectorXd y, 
                                Eigen::MatrixXd X, 
                                Eigen::MatrixXd Z,
                                Eigen::VectorXd beta, 
                                Eigen::VectorXd gamma, 
                                double sig2) {
  flexEL::MeanRegModel MR(y, X, Z);
  MatrixXd G = MatrixXd::Zero(MR.get_n_eqs(), MR.get_n_obs()); // G in C++ is n_eqs * n_obs
  MR.EvalG(G, beta, gamma, sig2);
  return(G);
}
