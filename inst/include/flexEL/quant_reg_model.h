/**
 * @file quant_reg_model.h
 * 
 * @brief A base class for template classes InnerEL and InnerELC.
 */

#ifndef QUANT_REG_MODEL_H
#define QUANT_REG_MODEL_H

#include "ind_smooth.h"

using namespace Eigen;

/* --------------------------------------------------------------------------- */

/**
* @brief EL namespace
* 
* Wrap the exported library components into a namespace called \b el to avoid potential naming conflicts with other libraries or user-defined headers.
*/
namespace flexEL {
  
  /**
  * @file       quant_reg_model.h
  * 
  * @class      QuantRegModel
  * 
  * @brief      A base class for template classes InnerEL and InnerELC to calculate estimating equations for quantile regressions.
  */
  class QuantRegModel {
  private:
    
    RowVectorXd yXb_;
    RowVectorXd eZg_;
    RowVectorXd yXbeZg_;
    RowVectorXd yXbeZg2_;
    MatrixXd tG_;
    int n_bet_, n_gam_, n_qts_; // nQts is the number of quantile levels
    VectorXd y_;
    MatrixXd X_;
    MatrixXd Z_;
    double *tau; // an array of quantile levels
    
    double phi_tau(double u, double tau); 
    double phi_tau_smooth(double u, double tau, double s); 
    
  protected:
    
    int n_obs_; /**< number of observations (number of columns of G) */
    int n_eqs_; /**< number of estimating equations (number of rows of G) */
  
  public:
    
    // constructors
    QuantRegModel();
    QuantRegModel(int n_obs, int n_eqs);
    QuantRegModel(const Ref<const VectorXd>& y,
                  const Ref<const MatrixXd>& X,
                  void* params);
    QuantRegModel(const Ref<const VectorXd>& y,
                  const Ref<const MatrixXd>& X,
                  const Ref<const MatrixXd>& Z,
                  void* params);
    
    // set data functions
    void set_data(const Ref<const VectorXd>& y, 
                  const Ref<const MatrixXd>& X,
                  void* params); // set data with default ctor
    void set_data(const Ref<const VectorXd>& y, 
                  const Ref<const MatrixXd>& X,
                  const Ref<const MatrixXd>& Z,
                  void* params); // set data with default ctor
    
    // evaluate G matrix
    void EvalG(Ref<MatrixXd> G, const Ref<const MatrixXd>& Beta);
    void EvalG(Ref<MatrixXd> G, 
               const Ref<const VectorXd>& beta, 
               const Ref<const VectorXd>& gamma,
               const double& sig2, // sig2 should be a scalar
               const Ref<const VectorXd>& nu);
    void EvalGSmooth(Ref<MatrixXd> G, 
                     const Ref<const MatrixXd>& Beta, 
                     const double s);
    void EvalGSmooth(Ref<MatrixXd> G,
                     const Ref<const VectorXd>& beta,
                     const Ref<const VectorXd>& gamma,
                     const double& sig2, // sig2 should be a scalar
                     const Ref<const VectorXd>& nu,
                     const double s);
    
    // get functions
    int get_n_obs();
    int get_n_eqs();
    
  };
}

/* --------------------------------------------------------------------------- */

// private functions

// 1st derivative of rho_tau
/**
* @brief First derivative of the check function.
* 
* @param[in] u     Argument of check function.
* @param[in] tau   Quantile level (0 < tau < 1).
*/
inline double flexEL::QuantRegModel::phi_tau(double u, double tau) {
  return((u <= 0) - tau);
}

// 1st derivative of rho_tau_smooth
/**
* @brief First derivative of the smoothed check function.
* 
* @param[in] u     Argument of check function.
* @param[in] tau   Quantile level  (0 < tau < 1).
* @param[in] s     Smoothing parameter (s > 0).
*/
inline double flexEL::QuantRegModel::phi_tau_smooth(double u, double tau, double s) {
  return(tau-ind_smooth(u,s)-u*ind1_smooth(u,s));
}

/* --------------------------------------------------------------------------- */

// public functions

// ctors
/**
* @brief Default constructor for QuantRegModel.
*/
inline flexEL::QuantRegModel::QuantRegModel(){}

/**
* @brief Constructor for QuantRegModel with dimensions as inputs.
* 
* @param[in] n_obs    Number of observations.
* @param[in] n_eqs    Number of estimating equations.
*/
inline flexEL::QuantRegModel::QuantRegModel(int n_obs, int n_eqs) {
  n_obs_ = n_obs;
  n_eqs_ = n_eqs; // X gets passed as nBet x n_obs matrix
}

/**
 * @brief Constructor for QuantRegModel with data as inputs (location model).
 * 
 * @param[in] y      Responses of length <code>n_obs</code>.
 * @param[in] X      Covariate matrix of dimension <code>nBet</code> x <code>n_obs</code>.
 * @param[in] params   An array of quantile levels.
 */
inline flexEL::QuantRegModel::QuantRegModel(const Ref<const VectorXd>& y,
                                            const Ref<const MatrixXd>& X,
                                            void* params) {
  y_ = y;
  X_ = X;
  double *tauArr = (double*)(params);
  n_qts_ = tauArr[0]; // first entry must be the number of quantile levels
  tau = &(tauArr[1]); // rest are the quantile levels
  
  n_obs_ = y.size();
  n_bet_ = X.rows(); // X gets passed as nBet x n_obs matrix
  n_gam_ = 0; 
  n_eqs_ = n_bet_*n_qts_; // Total number of equations
}

/**
 * @brief Constructor for QuantRegModel with data as inputs (location-scale model).
 * 
 * @param[in] y      Responses of length <code>n_obs</code>.
 * @param[in] X      Covariate matrix of dimension <code>nBet</code> x <code>n_obs</code>.
 * @param[in] Z      Covariate matrix of dimension <code>nGam</code> x <code>n_obs</code>.
 * @param[in] params   An array of quantile levels.
 */
inline flexEL::QuantRegModel::QuantRegModel(const Ref<const VectorXd>& y,
                                            const Ref<const MatrixXd>& X,
                                            const Ref<const MatrixXd>& Z,
                                            void* params) {
  y_ = y;
  X_ = X;
  Z_ = Z;
  double *tauArr = (double*)(params);
  n_qts_ = tauArr[0]; // first entry must be the number of quantile levels
  tau = &(tauArr[1]); // rest are the quantile levels
  
  n_obs_ = y.size();
  n_bet_ = X.rows(); // X gets passed as nBet x n_obs matrix
  n_gam_ = Z.rows(); // Z gets passed as nGam x n_obs matrix
  n_eqs_ = n_bet_+n_gam_+1+n_qts_; 
}

// set_data (with default ctor)
/**
* @brief Set data for quantile regression location model.
* 
* @param[in] y      Responses of length <code>n_obs</code>.
* @param[in] X      Covariate matrix of dimension <code>nBet</code> x <code>n_obs</code>.
* @param[in] params   An array of quantile levels.
*/
inline void flexEL::QuantRegModel::set_data(const Ref<const VectorXd>& y,
                                            const Ref<const MatrixXd>& X,
                                            void* params) {
  y_ = y;
  X_ = X;
  double *tauArr = (double*)(params);
  n_qts_ = tauArr[0]; // first entry must be the number of quantile levels
  tau = &(tauArr[1]); // rest are the quantile levels
  
  n_obs_ = y.size();
  n_bet_ = X.rows(); // X gets passed as nBet x n_obs matrix
  n_gam_ = 0; 
  n_eqs_ = n_bet_*n_qts_; // Total number of equations
}

// set_data (location-scale model) (with default ctor)
/**
* @brief Set data for quantile regression location-scale model.
* 
* @param[in] y      Responses of length <code>n_obs</code>.
* @param[in] X      Covariate matrix of dimension <code>nBet</code> x <code>n_obs</code>.
* @param[in] Z      Covariate matrix of dimension <code>nGam</code> x <code>n_obs</code>.
* @param[in] params   An array of quantile levels.
*/
inline void flexEL::QuantRegModel::set_data(const Ref<const VectorXd>& y,
                                            const Ref<const MatrixXd>& X,
                                            const Ref<const MatrixXd>& Z,
                                            void* params) {
  y_ = y;
  X_ = X;
  Z_ = Z;
  double *tauArr = (double*)(params);
  n_qts_ = tauArr[0]; // first entry must be the number of quantile levels
  tau = &(tauArr[1]); // rest are the quantile levels
  
  n_obs_ = y.size();
  n_bet_ = X.rows(); // X gets passed as nBet x n_obs matrix
  n_gam_ = Z.rows(); // Z gets passed as nGam x n_obs matrix
  n_eqs_ = n_bet_+n_gam_+1+n_qts_; 
}

// multiple quantile case (location model)
/**
* @brief Evaluate G matrix for quantile regression location model.
* 
* @param[in] G        Matrix of dimension <code>n_eqs_ x n_obs_</code> where the calculated result is saved.
* @param[in] Beta     Matrix of dimension <code>n_eqs_ x n_qts_</code>, each column is a coefficient vector of length <code>n_eqs_</code> in linear location function.
*/
inline void flexEL::QuantRegModel::EvalG(Ref<MatrixXd> G, const Ref<const MatrixXd>& Beta) {
  // stack multiple "Gs" for each quantile level
  for(int jj=0; jj<n_qts_; jj++) {
    for(int ii=0; ii<y_.size(); ii++) {
      G.block(jj*n_bet_,ii,n_bet_,1) = phi_tau(y_(ii)-X_.col(ii).transpose()*Beta.col(jj), tau[jj])*X_.col(ii);
    }
  }
}

// multiple quantile case (location-scale model)
// Note: only nu's depend on the quantile level
/**
* @brief Evaluate G matrix for quantile regression location-scale model.
* 
* @param[in] G        Matrix of dimension <code>n_eqs_ x n_obs_</code> where the calculated result is saved.
* @param[in] beta     Coefficient vector of length <code>nBet</code> in linear location function.
* @param[in] gamma    Coefficient vector of length <code>nGam</code> in exponential scale function.
* @param[in] sig2     Scale parameter in scale function.
* @param[in] nu       Quantile parameters for each quantile level.
*/
inline void flexEL::QuantRegModel::EvalG(Ref<MatrixXd> G, 
                                         const Ref<const VectorXd>& beta,
                                         const Ref<const VectorXd>& gamma,
                                         const double& sig2,
                                         const Ref<const VectorXd>& nu) {
  // TODO: warning if negative sig2?
  // if (sig2 < 0) std::cout << "qrls.EvalG: negative variance." << std::endl;
  eZg_.array() = (-gamma.transpose()*Z_).array().exp();
  yXbeZg_.array() = (y_.transpose()-beta.transpose()*X_).array() * eZg_.array();
  yXbeZg2_.array() = yXbeZg_.array()*yXbeZg_.array();
  
  tG_ = MatrixXd::Zero(n_obs_,n_eqs_); // NEW: DEC 25
  // 1st deriv w.r.t beta
  tG_.block(0,0,n_obs_,n_bet_) = X_.transpose();
  tG_.block(0,0,n_obs_,n_bet_).array().colwise() *= yXbeZg_.transpose().array() * eZg_.transpose().array();
  // 1st deriv w.r.t gamma
  tG_.block(0,n_bet_,n_obs_,n_gam_) = Z_.transpose();
  tG_.block(0,n_bet_,n_obs_,n_gam_).array().colwise() *= (1.0-yXbeZg2_.transpose().array()/sig2); // NEW (DEC 6)
  // variance param
  tG_.block(0,n_bet_+n_gam_,n_obs_,1).array() = 1/sig2*yXbeZg2_.transpose().array()-1;
  // quantile param(s)
  for (int ii=0; ii<n_obs_; ii++) {
    for (int jj=0; jj<n_qts_; jj++) {
      tG_.block(ii,n_bet_+n_gam_+1+jj,1,1).array() = phi_tau(yXbeZg_(ii)/sqrt(sig2)-nu(jj), tau[jj]);
    }
  }
  G = tG_.transpose();
}

// multiple quantile case (location model)
/**
 * @brief Evaluate G matrix for quantile regression location model.
 * 
 * @param[in] G        Matrix of dimension <code>n_eqs_ x n_obs_</code> where the calculated result is saved.
 * @param[in] Beta     Matrix of dimension <code>n_eqs_ x n_qts_</code>, each column is a coefficient vector of length <code>n_eqs_</code> in linear location function.
 */
inline void flexEL::QuantRegModel::EvalGSmooth(Ref<MatrixXd> G, 
                                               const Ref<const MatrixXd>& Beta,
                                               const double s) {
  // stack multiple "Gs" for each quantile level
  // std::cout << "n_bet_ = " << n_bet_ << std::endl;
  for(int jj=0; jj<n_qts_; jj++) {
    // std::cout << "** jj = " << jj << std::endl;
    for(int ii=0; ii<y_.size(); ii++) {
      // std::cout << "y_(ii) = " << y_(ii) << std::endl;
      // std::cout << "X_.col(ii) = " << X_.col(ii).transpose() << std::endl;
      // std::cout << "Beta.col(jj) = " << Beta.col(jj).transpose() << std::endl;
      // std::cout << "s = " << s << std::endl;
      // std::cout << "phi_tau_smooth(y_(ii)-X_.col(ii).transpose()*Beta.col(jj), tau[jj], s) = " << phi_tau_smooth(y_(ii)-X_.col(ii).transpose()*Beta.col(jj), tau[jj], s) << std::endl;
      G.block(jj*n_bet_,ii,n_bet_,1) = phi_tau_smooth(y_(ii)-X_.col(ii).transpose()*Beta.col(jj), tau[jj], s)*X_.col(ii);
    }
  }
}

/**
* @brief Evaluate G matrix for quantile regression location-scale model with smoothing.
* 
* @param[in] G        Matrix of dimension <code>n_eqs_ x n_obs_</code> where the calculated result is saved.
* @param[in] beta     Coefficient vector of length <code>nBet</code> in linear location function.
* @param[in] gamma    Coefficient vector of length <code>nGam</code> in exponential scale function.
* @param[in] sig2     Scale parameter in scale function.
* @param[in] nu       Quantile parameters for each quantile level.
* @param[in] s        Smoothing parameter (s > 0).
*/
inline void flexEL::QuantRegModel::EvalGSmooth(Ref<MatrixXd> G, 
                                               const Ref<const VectorXd>& beta,
                                               const Ref<const VectorXd>& gamma,
                                               const double& sig2,
                                               const Ref<const VectorXd>& nu,
                                               const double s) {
  // TODO: warning if negative sig2?
  // if (sig2 < 0) std::cout << "qrls.EvalG: negative variance." << std::endl;
  eZg_.array() = (-gamma.transpose()*Z_).array().exp();
  yXbeZg_.array() = (y_.transpose()-beta.transpose()*X_).array() * eZg_.array();
  yXbeZg2_.array() = yXbeZg_.array()*yXbeZg_.array();
  
  tG_ = MatrixXd::Zero(n_obs_,n_eqs_);
  // 1st deriv w.r.t beta
  tG_.block(0,0,n_obs_,n_bet_) = X_.transpose();
  tG_.block(0,0,n_obs_,n_bet_).array().colwise() *= yXbeZg_.transpose().array() * eZg_.transpose().array();
  // 1st deriv w.r.t gamma
  tG_.block(0,n_bet_,n_obs_,n_gam_) = Z_.transpose();
  tG_.block(0,n_bet_,n_obs_,n_gam_).array().colwise() *= (1.0-yXbeZg2_.transpose().array()/sig2); // NEW (DEC 6)
  // variance param
  tG_.block(0,n_bet_+n_gam_,n_obs_,1).array() = 1/sig2*yXbeZg2_.transpose().array()-1;
  // quantile param(s)
  for (int ii=0; ii<n_obs_; ii++) {
    for (int jj=0; jj<n_qts_; jj++) {
      tG_.block(ii,n_bet_+n_gam_+1+jj,1,1).array() = phi_tau_smooth(yXbeZg_(ii)/sqrt(sig2)-nu(jj),tau[jj],s);
    }
  }
  G = tG_.transpose();
}

/**
 * @brief Return the number of observations.
 */
inline int flexEL::QuantRegModel::get_n_obs() {
  return n_obs_;
}

/**
 * @brief Return the number of estimating equations.
 */
inline int flexEL::QuantRegModel::get_n_eqs() {
  return n_eqs_;
}

#endif
