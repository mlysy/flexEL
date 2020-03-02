/**
 * @file QuantRegModel.h
 * 
 * @brief A base class for template classes InnerEL and InnerELC.
 */

#ifndef QUANTREGMODEL_h
#define QUANTREGMODEL_h

#include "IndSmooth.h"

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
  int nBet_, nGam_, nQts_; // nQts is the number of quantile levels
  VectorXd y_;
  MatrixXd X_;
  MatrixXd Z_;
  double *tau; // an array of quantile levels
  
  double phi_tau(double u, double tau); 
  double phi_tau_smooth(double u, double tau, double s); 
  
protected:
  
  int nObs_; /**< number of observations (number of columns of G) */
  int nEqs_; /**< number of estimating equations (number of rows of G) */

public:
  
  // constructors
  QuantRegModel();
  QuantRegModel(int nObs, int nEqs);
  QuantRegModel(const Ref<const VectorXd>& y,
                const Ref<const MatrixXd>& X,
                void* params);
  QuantRegModel(const Ref<const VectorXd>& y,
                const Ref<const MatrixXd>& X,
                const Ref<const MatrixXd>& Z,
                void* params);
  
  // set data functions
  void setData(const Ref<const VectorXd>& y, 
               const Ref<const MatrixXd>& X,
               void* params); // set data with default ctor
  void setData(const Ref<const VectorXd>& y, 
               const Ref<const MatrixXd>& X,
               const Ref<const MatrixXd>& Z,
               void* params); // set data with default ctor
  
  // evaluate G matrix
  void evalG(Ref<MatrixXd> G, const Ref<const MatrixXd>& Beta);
  void evalG(Ref<MatrixXd> G, 
             const Ref<const VectorXd>& beta, 
             const Ref<const VectorXd>& gamma,
             const double& sig2, // sig2 should be a scalar
             const Ref<const VectorXd>& Nu);
  void evalGSmooth(Ref<MatrixXd> G,
                   const Ref<const VectorXd>& beta,
                   const Ref<const VectorXd>& gamma,
                   const double& sig2, // sig2 should be a scalar
                   const Ref<const VectorXd>& Nu,
                   const double s);
  
  // get functions
  int getnObs();
  int getnEqs();
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
 * @param nObs    Number of observations.
 * @param nEqs    Number of estimating equations.
 */
inline QuantRegModel::QuantRegModel(int nObs, int nEqs) {
  nObs_ = nObs;
  nEqs_ = nEqs; // X gets passed as nBet x nObs matrix
}

inline QuantRegModel::QuantRegModel(const Ref<const VectorXd>& y,
                                    const Ref<const MatrixXd>& X,
                                    void* params) {
  y_ = y;
  X_ = X;
  double *tauArr = (double*)(params);
  nQts_ = tauArr[0]; // first entry must be the number of quantile levels
  tau = &(tauArr[1]); // rest are the quantile levels
  
  nObs_ = y.size();
  nBet_ = X.rows(); // X gets passed as nBet x nObs matrix
  nGam_ = 0; 
  nEqs_ = nBet_*nQts_; // Total number of equations
}

inline QuantRegModel::QuantRegModel(const Ref<const VectorXd>& y,
                                    const Ref<const MatrixXd>& X,
                                    const Ref<const MatrixXd>& Z,
                                    void* params) {
  y_ = y;
  X_ = X;
  Z_ = Z;
  double *tauArr = (double*)(params);
  nQts_ = tauArr[0]; // first entry must be the number of quantile levels
  tau = &(tauArr[1]); // rest are the quantile levels
  
  nObs_ = y.size();
  nBet_ = X.rows(); // X gets passed as nBet x nObs matrix
  nGam_ = Z.rows(); // Z gets passed as nGam x nObs matrix
  nEqs_ = nBet_+nGam_+1+nQts_; 
}

// setData (with default ctor)
/**
 * @brief Set data for quantile regression location model.
 * 
 * @param y      Responses of length <code>nObs</code>.
 * @param X      Covariate matrix of dimension <code>nBet</code> x <code>nObs</code>.
 * @param params   An array of quantile levels.
 */
inline void QuantRegModel::setData(const Ref<const VectorXd>& y,
                                   const Ref<const MatrixXd>& X,
                                   void* params) {
  y_ = y;
  X_ = X;
  double *tauArr = (double*)(params);
  nQts_ = tauArr[0]; // first entry must be the number of quantile levels
  tau = &(tauArr[1]); // rest are the quantile levels
  
  nObs_ = y.size();
  nBet_ = X.rows(); // X gets passed as nBet x nObs matrix
  nGam_ = 0; 
  nEqs_ = nBet_*nQts_; // Total number of equations
}

// setData (location-scale model) (with default ctor)
/**
 * @brief Set data for quantile regression location-scale model.
 * 
 * @param y      Responses of length <code>nObs</code>.
 * @param X      Covariate matrix of dimension <code>nBet</code> x <code>nObs</code>.
 * @param Z      Covariate matrix of dimension <code>nGam</code> x <code>nObs</code>.
 * @param params   An array of quantile levels.
 */
inline void QuantRegModel::setData(const Ref<const VectorXd>& y,
                                   const Ref<const MatrixXd>& X,
                                   const Ref<const MatrixXd>& Z,
                                   void* params) {
  y_ = y;
  X_ = X;
  Z_ = Z;
  double *tauArr = (double*)(params);
  nQts_ = tauArr[0]; // first entry must be the number of quantile levels
  tau = &(tauArr[1]); // rest are the quantile levels
  
  nObs_ = y.size();
  nBet_ = X.rows(); // X gets passed as nBet x nObs matrix
  nGam_ = Z.rows(); // Z gets passed as nGam x nObs matrix
  nEqs_ = nBet_+nGam_+1+nQts_; 
}

// multiple quantile case (location model)
/**
 * @brief Evaluate G matrix for quantile regression location model.
 * 
 * @param Beta     Matrix of dimension <code>nEqs_ x nQts_</code>, each column is a coefficient vector of length <code>nEqs_</code> in linear location function.
 */
inline void QuantRegModel::evalG(Ref<MatrixXd> G, const Ref<const MatrixXd>& Beta) {
  // stack multiple "Gs" for each quantile level
  for(int jj=0; jj<nQts_; jj++) {
    for(int ii=0; ii<y_.size(); ii++) {
      G.block(jj*nBet_,ii,nBet_,1) = phi_tau(y_(ii)-X_.col(ii).transpose()*Beta.col(jj), tau[jj])*X_.col(ii);
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
inline void QuantRegModel::evalG(Ref<MatrixXd> G, 
                                 const Ref<const VectorXd>& beta,
                                 const Ref<const VectorXd>& gamma,
                                 const double& sig2,
                                 const Ref<const VectorXd>& Nu) {
  // TODO: warning if negative sig2?
  if (sig2 < 0) std::cout << "qrls.evalG: negative variance." << std::endl;
  eZg_.array() = (-gamma.transpose()*Z_).array().exp();
  yXbeZg_.array() = (y_.transpose()-beta.transpose()*X_).array() * eZg_.array();
  yXbeZg2_.array() = yXbeZg_.array()*yXbeZg_.array();
  
  tG_ = MatrixXd::Zero(nObs_,nEqs_); // NEW: DEC 25
  // 1st deriv w.r.t beta
  tG_.block(0,0,nObs_,nBet_) = X_.transpose();
  tG_.block(0,0,nObs_,nBet_).array().colwise() *= yXbeZg_.transpose().array() * eZg_.transpose().array();
  // 1st deriv w.r.t gamma
  tG_.block(0,nBet_,nObs_,nGam_) = Z_.transpose();
  tG_.block(0,nBet_,nObs_,nGam_).array().colwise() *= (1.0-yXbeZg2_.transpose().array()/sig2); // NEW (DEC 6)
  // variance param
  tG_.block(0,nBet_+nGam_,nObs_,1).array() = 1/sig2*yXbeZg2_.transpose().array()-1;
  // quantile param(s)
  for (int ii=0; ii<nObs_; ii++) {
    for (int jj=0; jj<nQts_; jj++) {
      tG_.block(ii,nBet_+nGam_+1+jj,1,1).array() = phi_tau(yXbeZg_(ii)/sqrt(sig2)-Nu(jj), tau[jj]);
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
inline void QuantRegModel::evalGSmooth(Ref<MatrixXd> G, 
                                       const Ref<const VectorXd>& beta,
                                       const Ref<const VectorXd>& gamma,
                                       const double& sig2,
                                       const Ref<const VectorXd>& Nu,
                                       const double s) {
  // TODO: warning if negative sig2?
  if (sig2 < 0) std::cout << "qrls.evalG: negative variance." << std::endl;
  eZg_.array() = (-gamma.transpose()*Z_).array().exp();
  yXbeZg_.array() = (y_.transpose()-beta.transpose()*X_).array() * eZg_.array();
  yXbeZg2_.array() = yXbeZg_.array()*yXbeZg_.array();
  
  tG_ = MatrixXd::Zero(nObs_,nEqs_);
  // 1st deriv w.r.t beta
  tG_.block(0,0,nObs_,nBet_) = X_.transpose();
  tG_.block(0,0,nObs_,nBet_).array().colwise() *= yXbeZg_.transpose().array() * eZg_.transpose().array();
  // 1st deriv w.r.t gamma
  tG_.block(0,nBet_,nObs_,nGam_) = Z_.transpose();
  tG_.block(0,nBet_,nObs_,nGam_).array().colwise() *= (1.0-yXbeZg2_.transpose().array()/sig2); // NEW (DEC 6)
  // variance param
  tG_.block(0,nBet_+nGam_,nObs_,1).array() = 1/sig2*yXbeZg2_.transpose().array()-1;
  // quantile param(s)
  for (int ii=0; ii<nObs_; ii++) {
    for (int jj=0; jj<nQts_; jj++) {
      tG_.block(ii,nBet_+nGam_+1+jj,1,1).array() = phi_tau_smooth(yXbeZg_(ii)/sqrt(sig2)-Nu(jj),tau[jj],s);
    }
  }
  G = tG_.transpose();
}

inline int QuantRegModel::getnObs() {
  return nObs_;
}

inline int QuantRegModel::getnEqs() {
  return nEqs_;
}

#endif