/**
 * @file QuantRegModel.h
 * 
 * @brief A base class for template classes InnerEL and InnerELC.
 */

#ifndef QUANTREGMODEL_h
#define QUANTREGMODEL_h

// #include <math.h>
// #include <Rmath.h>

#include "IndSmooth.h"

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
  
  /**
   * @brief First derivative of the check function.
   * 
   * @param u     Argument of check function.
   * @param tau   Quantile level (0 < tau < 1).
   */
  double phi_tau(double u, double tau); 
  
protected:
  
  int nObs_, nEqs_;
  int nBet_, nGam_, nQts_; /**< nQts is the number of quantile levels */ 
  VectorXd y_;
  MatrixXd X_;
  MatrixXd Z_;
  MatrixXd G_;
  double *tau; /**< an array of quantile levels */ 
  
public:
  
  /**
   * @brief Default constructor for QuantRegModel.
   */
  QuantRegModel();
  
  /**
   * @brief Constructor for QuantRegModel with dimensions as inputs.
   * 
   * @param nObs    Number of observations.
   * @param nEqs    Number of estimating equations.
   */
  QuantRegModel(int nObs, int nEqs);
  
  /**
   * @brief Set data for quantile regression location model.
   * 
   * @param y      Responses of length <code>nObs</code>.
   * @param X      Covariate matrix of dimension <code>nBet</code> x <code>nObs</code>.
   * @param params   An array of quantile levels.
   */
  void setData(const Ref<const VectorXd>& y, 
               const Ref<const MatrixXd>& X,
               void* params); // set data with default ctor
  
  /**
   * @brief Set data for quantile regression location-scale model.
   * 
   * @param y      Responses of length <code>nObs</code>.
   * @param X      Covariate matrix of dimension <code>nBet</code> x <code>nObs</code>.
   * @param Z      Covariate matrix of dimension <code>nGam</code> x <code>nObs</code>.
   * @param params   An array of quantile levels.
   */
  void setData(const Ref<const VectorXd>& y, 
               const Ref<const MatrixXd>& X,
               const Ref<const MatrixXd>& Z,
               void* params); // set data with default ctor
  
  // evaluate G matrix
  /**
   * @brief Evaluate G matrix for quantile regression location model.
   * 
   * @param beta     Coefficient vector of length <code>nObs</code> in linear location function.
   */
  void evalG(const Ref<const MatrixXd>& Beta);
  
  /**
   * @brief Evaluate G matrix for quantile regression location-scale model.
   * 
   * @param beta     Coefficient vector of length <code>nBet</code> in linear location function.
   * @param gamma    Coefficient vector of length <code>nGam</code> in exponential scale function.
   * @param sig2     Scale parameter in scale function.
   * @param Nu       Quantile parameters for each quantile level.
   */
  void evalG(const Ref<const VectorXd>& beta, 
             const Ref<const VectorXd>& gamma,
             const double& sig2, // sig2 should be a scalar
             const Ref<const VectorXd>& Nu);
  
  // evaluate G matrix (smoothed version)
  /**
   * @brief First derivative of the smoothed check function.
   * 
   * @param u     Argument of check function.
   * @param tau   Quantile level  (0 < tau < 1).
   * @param s     Smoothing parameter (s > 0).
   */
  double phi_tau_smooth(double u, double tau, double s); 
  
  /**
   * @brief Evaluate G matrix for smoothed quantile regression location-scale model.
   * 
   * @param beta     Coefficient vector of length <code>nObs</code> in linear location function.
   */
  void evalGSmooth(const Ref<const VectorXd>& beta,
                   const Ref<const VectorXd>& gamma,
                   const double& sig2, // sig2 should be a scalar
                   const Ref<const VectorXd>& Nu,
                   const double s);
};

// default ctor 
inline QuantRegModel::QuantRegModel(){}

// ctor
inline QuantRegModel::QuantRegModel(int nObs, int nEqs) {
  nObs_ = nObs;
  nEqs_ = nEqs; // X gets passed as nBet x nObs matrix
  G_ = MatrixXd::Zero(nEqs_,nObs_);
}

// setData (with default ctor)
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
  
  G_ = MatrixXd::Zero(nEqs_,nObs_);
  tG_ = MatrixXd::Zero(nObs_,nEqs_);
}

// setData (location-scale model) (with default ctor)
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
  // nEqs = nQts*(nBet+nGam+2);
  // nEqs = nQts*(nBet+nGam);
  G_ = MatrixXd::Zero(nEqs_,nObs_);
  tG_ = MatrixXd::Zero(nObs_,nEqs_);
  eZg_ = RowVectorXd::Zero(nObs_);
  yXbeZg_ = RowVectorXd::Zero(nObs_);
  yXbeZg2_ = RowVectorXd::Zero(nObs_);
}

// 1st derivative of rho_tau
inline double QuantRegModel::phi_tau(double u, double tau) {
    return((u <= 0) - tau);
}

// 1st derivative of rho_tau_smooth
inline double QuantRegModel::phi_tau_smooth(double u, double tau, double s) {
  return(tau-ind_smooth(u,s)-u*ind1_smooth(u,s));
}

/*
// single quantile case (location model)
inline void QuantRegModel::evalG(const Ref<const VectorXd>& beta) {
  for(int ii=0; ii<y.size(); ii++) {
      this->G.col(ii) = phi_tau(y(ii)-X.col(ii).transpose()*beta, tau[0])*X.col(ii);
  }
}
*/

// multiple quantile case (location model)
inline void QuantRegModel::evalG(const Ref<const MatrixXd>& Beta) {
  // stack multiple "Gs" for each quantile level
  for(int jj=0; jj<nQts_; jj++) {
    for(int ii=0; ii<y_.size(); ii++) {
      this->G_.block(jj*nBet_,ii,nBet_,1) = phi_tau(y_(ii)-X_.col(ii).transpose()*Beta.col(jj), tau[jj])*X_.col(ii);
    }
  }
}

// multiple quantile case (location-scale model)
// Note: only nu's depend on the quantile level
inline void QuantRegModel::evalG(const Ref<const VectorXd>& beta,
                                 const Ref<const VectorXd>& gamma,
                                 const double& sig2,
                                 const Ref<const VectorXd>& Nu) {
  // std::cout << "nQts = " << nQts << std::endl;
  // TODO: warning if negative sig2?
  if (sig2 < 0) std::cout << "qrls.evalG: negative variance." << std::endl;
  eZg_.array() = (-gamma.transpose()*Z_).array().exp();
  yXbeZg_.array() = (y_.transpose()-beta.transpose()*X_).array() * eZg_.array();
  yXbeZg2_.array() = yXbeZg_.array()*yXbeZg_.array();
  // 1st deriv w.r.t beta
  tG_.block(0,0,nObs_,nBet_) = X_.transpose();
  tG_.block(0,0,nObs_,nBet_).array().colwise() *= yXbeZg_.transpose().array() * eZg_.transpose().array();
  // 1st deriv w.r.t gamma
  tG_.block(0,nBet_,nObs_,nGam_) = Z_.transpose();
  // tG.block(0,nBet,nObs,nGam).array().colwise() *= yXbeZg2.transpose().array();
  // tG.block(0,nBet,nObs,nGam).array().colwise() *= (1.0-yXbeZg2.transpose().array());
  tG_.block(0,nBet_,nObs_,nGam_).array().colwise() *= (1.0-yXbeZg2_.transpose().array()/sig2); // NEW (DEC 6)
  // variance param
  tG_.block(0,nBet_+nGam_,nObs_,1).array() = 1/sig2*yXbeZg2_.transpose().array()-1;
  // std::cout << "tG = \n" << tG << std::endl;
  // quantile param(s)
  for (int ii=0; ii<nObs_; ii++) {
    for (int jj=0; jj<nQts_; jj++) {
      tG_.block(ii,nBet_+nGam_+1+jj,1,1).array() = phi_tau(yXbeZg_(ii)/sqrt(sig2)-Nu(jj), tau[jj]);
    }
  }
  G_ = tG_.transpose();
}


inline void QuantRegModel::evalGSmooth(const Ref<const VectorXd>& beta,
                                       const Ref<const VectorXd>& gamma,
                                       const double& sig2,
                                       const Ref<const VectorXd>& Nu,
                                       const double s) {
  // std::cout << "nBet = " << nBet << std::endl;
  // std::cout << "nGam = " << nGam << std::endl;
  // std::cout << "nrow = " << G.rows() << std::endl;
  // std::cout << "ncol = " << G.cols() << std::endl;
  // TODO: warning if negative sig2?
  if (sig2 < 0) std::cout << "qrls.evalG: negative variance." << std::endl;
  eZg_.array() = (-gamma.transpose()*Z_).array().exp();
  yXbeZg_.array() = (y_.transpose()-beta.transpose()*X_).array() * eZg_.array();
  yXbeZg2_.array() = yXbeZg_.array()*yXbeZg_.array();
  // 1st deriv w.r.t beta
  tG_.block(0,0,nObs_,nBet_) = X_.transpose();
  tG_.block(0,0,nObs_,nBet_).array().colwise() *= yXbeZg_.transpose().array() * eZg_.transpose().array();
  // 1st deriv w.r.t gamma
  tG_.block(0,nBet_,nObs_,nGam_) = Z_.transpose();
  // tG.block(0,nBet,nObs,nGam).array().colwise() *= yXbeZg2.transpose().array();
  // tG.block(0,nBet,nObs,nGam).array().colwise() *= (1.0-yXbeZg2.transpose().array());
  tG_.block(0,nBet_,nObs_,nGam_).array().colwise() *= (1.0-yXbeZg2_.transpose().array()/sig2); // NEW (DEC 6)
  // variance param
  tG_.block(0,nBet_+nGam_,nObs_,1).array() = 1/sig2*yXbeZg2_.transpose().array()-1;
  // std::cout << "tG = \n" << tG << std::endl;
  // quantile param(s)
  for (int ii=0; ii<nObs_; ii++) {
    for (int jj=0; jj<nQts_; jj++) {
      tG_.block(ii,nBet_+nGam_+1+jj,1,1).array() = phi_tau_smooth(yXbeZg_(ii)/sqrt(sig2)-Nu(jj),tau[jj],s);
    }
  }
  G_ = tG_.transpose();
}

#endif