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

/* ------ exported functions ------ */

// [[Rcpp::export(".MeanRegEvalG")]]
Eigen::MatrixXd MeanReg_evalG(Eigen::VectorXd y, Eigen::MatrixXd X, 
                              Eigen::VectorXd beta) {
  flexEL::MeanRegModel MR(y, X);
  MatrixXd G = MatrixXd::Zero(MR.get_n_eqs(), MR.get_n_obs()); // G in C++ is n_eqs * n_obs
  MR.EvalG(G, beta);
  return(G);
}

// [[Rcpp::export(".MeanRegLSEvalG")]]
Eigen::MatrixXd MeanRegLS_EvalG(Eigen::VectorXd y, 
                                Eigen::MatrixXd X, Eigen::MatrixXd Z,
                                Eigen::VectorXd beta, Eigen::VectorXd gamma, 
                                double sig2) {
  flexEL::MeanRegModel MR(y, X, Z);
  MatrixXd G = MatrixXd::Zero(MR.get_n_eqs(), MR.get_n_obs()); // G in C++ is n_eqs * n_obs
  MR.EvalG(G, beta, gamma, sig2);
  return(G);
}
