/**
 * @file quant_reg_model_exports.cpp
 * 
 * @brief Export quant_reg_model functions to R.
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include "inner_el.h"
#include "quant_reg_model.h"
#include "dvec_to_barr.h"

//[[Rcpp::depends("RcppEigen")]]

using namespace Rcpp;
using namespace Eigen;

/* ------ exported functions ------ */

// double tau
// [[Rcpp::export(".QuantRegEvalG")]]
Eigen::MatrixXd QuantRegEvalG(Eigen::VectorXd y, Eigen::MatrixXd X,
                               Eigen::VectorXd tauArr, Eigen::VectorXd beta) {
  flexEL::QuantRegModel QR(y,X,tauArr.data());
  flexEL::InnerEL EL(QR.get_n_obs(), QR.get_n_eqs());
  QR.EvalG(EL.get_ref_G(), beta);
  return(EL.get_G()); 
}

// [[Rcpp::export(".QuantRegLSEvalG")]]
Eigen::MatrixXd QuantRegLSEvalG(Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::MatrixXd Z, 
                                 Eigen::VectorXd tauArr, Eigen::VectorXd beta, 
                                 Eigen::VectorXd gamma, double sig2, Eigen::VectorXd nu) {
  flexEL::QuantRegModel QR(y,X,Z,tauArr.data());
  flexEL::InnerEL EL(QR.get_n_obs(), QR.get_n_eqs());
  QR.EvalG(EL.get_ref_G(),beta,gamma,sig2,nu);
  return(EL.get_G()); 
}

// [[Rcpp::export(".QuantRegLSEvalGSmooth")]]
Eigen::MatrixXd QuantRegLSEvalGSmooth(Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::MatrixXd Z, 
                                       Eigen::VectorXd tauArr, Eigen::VectorXd beta, 
                                       Eigen::VectorXd gamma, double sig2, Eigen::VectorXd nu, double s) {
  flexEL::QuantRegModel QR(y,X,Z,tauArr.data());
  flexEL::InnerEL EL(QR.get_n_obs(), QR.get_n_eqs());
  QR.EvalGSmooth(EL.get_ref_G(),beta,gamma,sig2,nu,s);
  return(EL.get_G()); 
}
