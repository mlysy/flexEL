/**
 * @file quant_reg_model_exports.cpp
 * 
 * @brief Export quant_reg_model functions to R.
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include "gen_el.h"
#include "quant_reg_model.h"

//[[Rcpp::depends("RcppEigen")]]

// using namespace Rcpp;
using namespace Eigen;

/* ------ exported functions ------ */

/**
 * @brief Calculate the G matrix for regular quantile regression model.
 * 
 * @param[in] y        Responses of length <code>n_obs</code>.
 * @param[in] X        Covariate matrix of dimension <code>nBet</code> x <code>n_obs</code>.
 * @param[in] tauArr   A numeric array of quantile levels.
 * @param[in] Beta     Coefficient vector of length <code>nBet</code> in linear location function.
 */
// [[Rcpp::export]]
Eigen::MatrixXd QuantReg_evalG(Eigen::VectorXd y, 
                               Eigen::MatrixXd X,
                               Eigen::VectorXd tauArr, 
                               Eigen::MatrixXd Beta) {
  flexEL::QuantRegModel QR(y, X, tauArr.data());
  MatrixXd G = MatrixXd::Zero(QR.get_n_eqs(), QR.get_n_obs()); // G in C++ is n_eqs * n_obs
  QR.eval_G(G, Beta);
  return(G);
}

/**
 * @brief Calculate the G matrix for location-scale quantile regression model.
 * 
 * @param[in] y      Responses of length <code>n_obs</code>.
 * @param[in] X      Covariate matrix of dimension <code>nBet</code> x <code>n_obs</code>.
 * @param[in] Z      Covariate matrix of dimension <code>nGam</code> x <code>n_obs</code>.
 * @param[in] tauArr   A numeric array of quantile levels.
 * @param[in] beta     Coefficient vector of length <code>nBet</code> in linear location function.
 * @param[in] gamma    Coefficient vector of length <code>nGam</code> in exponential scale function.
 * @param[in] sig2     Scale parameter in scale function.
 * @param[in] nu       Quantile parameters for each quantile level.
 */
// [[Rcpp::export]]
Eigen::MatrixXd QuantRegLS_evalG(Eigen::VectorXd y,
                                 Eigen::MatrixXd X, 
                                 Eigen::MatrixXd Z, 
                                 Eigen::VectorXd tauArr, 
                                 Eigen::VectorXd beta, 
                                 Eigen::VectorXd gamma, 
                                 double sig2, 
                                 Eigen::VectorXd nu) {
  flexEL::QuantRegModel QR(y, X, Z, tauArr.data());
  MatrixXd G = MatrixXd::Zero(QR.get_n_eqs(), QR.get_n_obs()); // G in C++ is n_eqs * n_obs
  QR.eval_G(G, beta, gamma, sig2, nu);
  return(G);
}

/**
 * @brief Calculate the G matrix for regular quantile regression model.
 * 
 * @param[in] y        Responses of length <code>n_obs</code>.
 * @param[in] X        Covariate matrix of dimension <code>nBet</code> x <code>n_obs</code>.
 * @param[in] tauArr   A numeric array of quantile levels.
 * @param[in] Beta     Coefficient vector of length <code>nBet</code> in linear location function.
 * @param[in] s        Tuning parameter for continuity correction.
 */
// [[Rcpp::export]]
Eigen::MatrixXd QuantReg_evalG_smooth(Eigen::VectorXd y, 
                                      Eigen::MatrixXd X,
                                      Eigen::VectorXd tauArr, 
                                      Eigen::MatrixXd Beta,
                                      double s) {
  flexEL::QuantRegModel QR(y, X, tauArr.data());
  // std::cout << "n_obs = " << QR.get_n_obs() << std::endl;
  // std::cout << "n_eqs = " << QR.get_n_eqs() << std::endl;
  MatrixXd G = MatrixXd::Zero(QR.get_n_eqs(), QR.get_n_obs()); // G in C++ is n_eqs * n_obs
  QR.eval_G_smooth(G, Beta, s);
  return(G);
}

/**
 * @brief Calculate the derivative of G matrix w.r.t Beta for regular quantile regression model with continuity correction.
 * 
 * @param[in] y        Responses of length <code>n_obs</code>.
 * @param[in] X        Covariate matrix of dimension <code>nBet</code> x <code>n_obs</code>.
 * @param[in] tauArr   A numeric array of quantile levels.
 * @param[in] Beta     Coefficient vector of length <code>nBet</code> in linear location function.
 * @param[in] s        Tuning parameter for continuity correction.
 */
// [[Rcpp::export]]
Eigen::MatrixXd QuantReg_dGdt_smooth(Eigen::VectorXd y, 
                                     Eigen::MatrixXd X,
                                     Eigen::VectorXd tauArr, 
                                     Eigen::MatrixXd Beta,
                                     double s) {
  flexEL::QuantRegModel QR(y, X, tauArr.data());
  // std::cout << "n_obs = " << QR.get_n_obs() << std::endl;
  // std::cout << "n_eqs = " << QR.get_n_eqs() << std::endl;
  MatrixXd dGdt = MatrixXd::Zero(QR.get_n_obs()*QR.get_n_bet(), QR.get_n_eqs());
  QR.eval_dGdt_smooth(dGdt, Beta, s);
  return(dGdt);
}

/**
 * @brief Calculate the G matrix for location-scale quantile regression model with continuity correction.
 * 
 * @param[in] y      Responses of length <code>n_obs</code>.
 * @param[in] X      Covariate matrix of dimension <code>nBet</code> x <code>n_obs</code>.
 * @param[in] Z      Covariate matrix of dimension <code>nGam</code> x <code>n_obs</code>.
 * @param[in] tauArr   A numeric array of quantile levels.
 * @param[in] beta     Coefficient vector of length <code>nBet</code> in linear location function.
 * @param[in] gamma    Coefficient vector of length <code>nGam</code> in exponential scale function.
 * @param[in] sig2     Scale parameter in scale function.
 * @param[in] nu       Quantile parameters for each quantile level.
 * @param[in] s        Smoothing parameter (s > 0).
 */
// [[Rcpp::export]]
Eigen::MatrixXd QuantRegLS_evalG_smooth(Eigen::VectorXd y,
                                        Eigen::MatrixXd X, 
                                        Eigen::MatrixXd Z, 
                                        Eigen::VectorXd tauArr, 
                                        Eigen::VectorXd beta, 
                                        Eigen::VectorXd gamma, 
                                        double sig2, 
                                        Eigen::VectorXd nu, 
                                        double s) {
  flexEL::QuantRegModel QR(y,X,Z,tauArr.data());
  MatrixXd G = MatrixXd::Zero(QR.get_n_eqs(), QR.get_n_obs()); // G in C++ is n_eqs * n_obs
  QR.eval_G_smooth(G, beta, gamma, sig2, nu, s);
  return(G);
}
