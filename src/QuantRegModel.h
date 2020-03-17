/**
 * @file QuantRegModel.h
 * 
 * @brief A base class for template classes InnerEL and InnerELC.
 */

#ifndef QUANTREGMODEL_h
#define QUANTREGMODEL_h

#include "ind_smooth.h"

/* --------------------------------------------------------------------------- */

/**
 * @file       QuantRegModel.h
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
             const Ref<const VectorXd>& Nu);
  void EvalGSmooth(Ref<MatrixXd> G,
                   const Ref<const VectorXd>& beta,
                   const Ref<const VectorXd>& gamma,
                   const double& sig2, // sig2 should be a scalar
                   const Ref<const VectorXd>& Nu,
                   const double s);
  
  // get functions
  int get_n_obs();
  int get_n_eqs();
};

/* --------------------------------------------------------------------------- */

// private functions

// 1st derivative of rho_tau
/**
* @brief First derivative of the check function.
* 
* @param u     Argument of check function.
* @param tau   Quantile level (0 < tau < 1).
*/
inline double QuantRegModel::phi_tau(double u, double tau) {
  return((u <= 0) - tau);
}

// 1st derivative of rho_tau_smooth
/**
* @brief First derivative of the smoothed check function.
* 
* @param u     Argument of check function.
* @param tau   Quantile level  (0 < tau < 1).
* @param s     Smoothing parameter (s > 0).
*/
inline double QuantRegModel::phi_tau_smooth(double u, double tau, double s) {
  return(tau-ind_smooth(u,s)-u*ind1_smooth(u,s));
}

/* --------------------------------------------------------------------------- */

// public functions

// ctors
/**
 * @brief Default constructor for QuantRegModel.
 */
inline QuantRegModel::QuantRegModel(){}

/**
 * @brief Constructor for QuantRegModel with dimensions as inputs.
 * 
 * @param n_obs    Number of observations.
 * @param n_eqs    Number of estimating equations.
 */
inline QuantRegModel::QuantRegModel(int n_obs, int n_eqs) {
  n_obs_ = n_obs;
  n_eqs_ = n_eqs; // X gets passed as nBet x n_obs matrix
}

inline QuantRegModel::QuantRegModel(const Ref<const VectorXd>& y,
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

inline QuantRegModel::QuantRegModel(const Ref<const VectorXd>& y,
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
 * @param y      Responses of length <code>n_obs</code>.
 * @param X      Covariate matrix of dimension <code>nBet</code> x <code>n_obs</code>.
 * @param params   An array of quantile levels.
 */
inline void QuantRegModel::set_data(const Ref<const VectorXd>& y,
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
 * @param y      Responses of length <code>n_obs</code>.
 * @param X      Covariate matrix of dimension <code>nBet</code> x <code>n_obs</code>.
 * @param Z      Covariate matrix of dimension <code>nGam</code> x <code>n_obs</code>.
 * @param params   An array of quantile levels.
 */
inline void QuantRegModel::set_data(const Ref<const VectorXd>& y,
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
 * @param Beta     Matrix of dimension <code>n_eqs_ x n_qts_</code>, each column is a coefficient vector of length <code>n_eqs_</code> in linear location function.
 */
inline void QuantRegModel::EvalG(Ref<MatrixXd> G, const Ref<const MatrixXd>& Beta) {
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
 * @param beta     Coefficient vector of length <code>nBet</code> in linear location function.
 * @param gamma    Coefficient vector of length <code>nGam</code> in exponential scale function.
 * @param sig2     Scale parameter in scale function.
 * @param Nu       Quantile parameters for each quantile level.
 */
inline void QuantRegModel::EvalG(Ref<MatrixXd> G, 
                                 const Ref<const VectorXd>& beta,
                                 const Ref<const VectorXd>& gamma,
                                 const double& sig2,
                                 const Ref<const VectorXd>& Nu) {
  // TODO: warning if negative sig2?
  if (sig2 < 0) std::cout << "qrls.EvalG: negative variance." << std::endl;
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
      tG_.block(ii,n_bet_+n_gam_+1+jj,1,1).array() = phi_tau(yXbeZg_(ii)/sqrt(sig2)-Nu(jj), tau[jj]);
    }
  }
  G = tG_.transpose();
}

/**
 * @brief Evaluate G matrix for quantile regression location-scale model with smoothing.
 * 
 * @param beta     Coefficient vector of length <code>nBet</code> in linear location function.
 * @param gamma    Coefficient vector of length <code>nGam</code> in exponential scale function.
 * @param sig2     Scale parameter in scale function.
 * @param Nu       Quantile parameters for each quantile level.
 * @param s        Smoothing parameter (s > 0).
 */
inline void QuantRegModel::EvalGSmooth(Ref<MatrixXd> G, 
                                       const Ref<const VectorXd>& beta,
                                       const Ref<const VectorXd>& gamma,
                                       const double& sig2,
                                       const Ref<const VectorXd>& Nu,
                                       const double s) {
  // TODO: warning if negative sig2?
  if (sig2 < 0) std::cout << "qrls.EvalG: negative variance." << std::endl;
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
      tG_.block(ii,n_bet_+n_gam_+1+jj,1,1).array() = phi_tau_smooth(yXbeZg_(ii)/sqrt(sig2)-Nu(jj),tau[jj],s);
    }
  }
  G = tG_.transpose();
}

inline int QuantRegModel::get_n_obs() {
  return n_obs_;
}

inline int QuantRegModel::get_n_eqs() {
  return n_eqs_;
}

#endif