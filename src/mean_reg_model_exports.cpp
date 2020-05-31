/**
 * @file mean_reg_model_exports.cpp
 * 
 * @brief Export mean_reg_model functions to R.
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include "inner_el.h"
#include "mean_reg_model.h"
#include "dvec_to_barr.h"

//[[Rcpp::depends("RcppEigen")]]

using namespace Rcpp;
using namespace Eigen;

/* ------ exported functions ------ */

// [[Rcpp::export(".MeanRegEvalG")]]
Eigen::MatrixXd MeanReg_evalG(Eigen::VectorXd y, Eigen::MatrixXd X, 
                              Eigen::VectorXd beta) {
  flexEL::MeanRegModel MR(y, X);
  flexEL::InnerEL EL(MR.get_n_obs(), MR.get_n_eqs());
  MR.EvalG(EL.get_ref_G(), beta);
  return(EL.get_G());
}

// [[Rcpp::export(".MeanRegLSEvalG")]]
Eigen::MatrixXd MeanRegLS_EvalG(Eigen::VectorXd y, 
                                Eigen::MatrixXd X, Eigen::MatrixXd Z,
                                Eigen::VectorXd beta, Eigen::VectorXd gamma, 
                                double sig2) {
  flexEL::MeanRegModel MR(y, X, Z);
  flexEL::InnerEL EL(MR.get_n_obs(), MR.get_n_eqs());
  MR.EvalG(EL.get_ref_G(), beta, gamma, sig2);
  return(EL.get_G());
}
