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
  
  RowVectorXd yXb;
  RowVectorXd eZg;
  RowVectorXd yXbeZg;
  RowVectorXd yXbeZg2;
  MatrixXd tG;
  
  /**
   * @brief First derivative of the check function.
   * 
   * @param u     Argument of check function.
   * @param tau   Quantile level (0 < tau < 1).
   */
  double phi_tau(double u, double tau); 
  
protected:
  
  int nObs, nEqs, nBet, nGam, nQts; /**< nQts is the number of quantile levels */ 
  VectorXd y;
  MatrixXd X;
  MatrixXd Z;
  MatrixXd G;
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
   * @param y      Responses of length \code{nObs}.
   * @param X      Covariate matrix of dimension \code{nBet} x \code{nObs}.
   * @parm params   An array of quantile levels.
   */
  void setData(const Ref<const VectorXd>& y, 
               const Ref<const MatrixXd>& X,
               void* params); // set data with default ctor
  
  /**
   * @brief Set data for quantile regression location-scale model.
   * 
   * @param y      Responses of length \code{nObs}.
   * @param X      Covariate matrix of dimension \code{nBet} x \code{nObs}.
   * @param Z      Covariate matrix of dimension \code{nGam} x \code{nObs}.
   * @parm params   An array of quantile levels.
   */
  void setData(const Ref<const VectorXd>& y, 
               const Ref<const MatrixXd>& X,
               const Ref<const MatrixXd>& Z,
               void* params); // set data with default ctor
  
  // evaluate G matrix
  /**
   * @brief Evaluate G matrix for quantile regression location model.
   * 
   * @param beta     Coefficient vector in linear location function.
   */
  void evalG(const Ref<const MatrixXd>& Beta);
  
  /**
   * @brief Evaluate G matrix for quantile regression location-scale model.
   * 
   * @param beta     Coefficient vector of length \code{nBet} in linear location function.
   * @param gamma    Coefficient vector of length \code{nGam} in exponential scale function.
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
   * @param beta     Coefficient vector in linear location function.
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
inline QuantRegModel::QuantRegModel(int _nObs, int _nEqs) {
  nObs = _nObs;
  nEqs = _nEqs; // X gets passed as nBet x nObs matrix
  G = MatrixXd::Zero(nEqs,nObs);
}

// setData (with default ctor)
inline void QuantRegModel::setData(const Ref<const VectorXd>& _y,
                                   const Ref<const MatrixXd>& _X,
                                   void* params) {
  y = _y;
  X = _X;
  double *tauArr = (double*)(params);
  nQts = tauArr[0]; // first entry must be the number of quantile levels
  tau = &(tauArr[1]); // rest are the quantile levels
  
  nObs = y.size();
  nBet = X.rows(); // X gets passed as nBet x nObs matrix
  nGam = 0; 
  nEqs = nBet*nQts; // Total number of equations
  
  G = MatrixXd::Zero(nEqs,nObs);
  tG = MatrixXd::Zero(nObs,nEqs);
}

// setData (location-scale model) (with default ctor)
inline void QuantRegModel::setData(const Ref<const VectorXd>& _y,
                                   const Ref<const MatrixXd>& _X,
                                   const Ref<const MatrixXd>& _Z,
                                   void* params) {
  y = _y;
  X = _X;
  Z = _Z;
  double *tauArr = (double*)(params);
  nQts = tauArr[0]; // first entry must be the number of quantile levels
  tau = &(tauArr[1]); // rest are the quantile levels
  
  nObs = y.size();
  nBet = X.rows(); // X gets passed as nBet x nObs matrix
  nGam = Z.rows(); // Z gets passed as nGam x nObs matrix
  nEqs = nBet+nGam+1+nQts; 
  // nEqs = nQts*(nBet+nGam+2);
  // nEqs = nQts*(nBet+nGam);
  G = MatrixXd::Zero(nEqs,nObs);
  tG = MatrixXd::Zero(nObs,nEqs);
  eZg = RowVectorXd::Zero(nObs);
  yXbeZg = RowVectorXd::Zero(nObs);
  yXbeZg2 = RowVectorXd::Zero(nObs);
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
  for(int jj=0; jj<nQts; jj++) {
    for(int ii=0; ii<y.size(); ii++) {
      this->G.block(jj*nBet,ii,nBet,1) = phi_tau(y(ii)-X.col(ii).transpose()*Beta.col(jj), tau[jj])*X.col(ii);
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
  eZg.array() = (-gamma.transpose()*Z).array().exp();
  yXbeZg.array() = (y.transpose()-beta.transpose()*X).array() * eZg.array();
  yXbeZg2.array() = yXbeZg.array()*yXbeZg.array();
  // 1st deriv w.r.t beta
  tG.block(0,0,nObs,nBet) = X.transpose();
  tG.block(0,0,nObs,nBet).array().colwise() *= yXbeZg.transpose().array() * eZg.transpose().array();
  // 1st deriv w.r.t gamma
  tG.block(0,nBet,nObs,nGam) = Z.transpose();
  // tG.block(0,nBet,nObs,nGam).array().colwise() *= yXbeZg2.transpose().array();
  tG.block(0,nBet,nObs,nGam).array().colwise() *= (1.0-yXbeZg2.transpose().array());
  // variance param
  tG.block(0,nBet+nGam,nObs,1).array() = 1/sig2*yXbeZg2.transpose().array()-1;
  // std::cout << "tG = \n" << tG << std::endl;
  // quantile param(s)
  for (int ii=0; ii<nObs; ii++) {
    for (int jj=0; jj<nQts; jj++) {
      tG.block(ii,nBet+nGam+1+jj,1,1).array() = phi_tau(yXbeZg(ii)/sqrt(sig2)-Nu(jj), tau[jj]);
    }
  }
  G = tG.transpose();
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
  eZg.array() = (-gamma.transpose()*Z).array().exp();
  yXbeZg.array() = (y.transpose()-beta.transpose()*X).array() * eZg.array();
  yXbeZg2.array() = yXbeZg.array()*yXbeZg.array();
  // 1st deriv w.r.t beta
  tG.block(0,0,nObs,nBet) = X.transpose();
  tG.block(0,0,nObs,nBet).array().colwise() *= yXbeZg.transpose().array() * eZg.transpose().array();
  // 1st deriv w.r.t gamma
  tG.block(0,nBet,nObs,nGam) = Z.transpose();
  // tG.block(0,nBet,nObs,nGam).array().colwise() *= yXbeZg2.transpose().array();
  tG.block(0,nBet,nObs,nGam).array().colwise() *= (1.0-yXbeZg2.transpose().array());
  // variance param
  tG.block(0,nBet+nGam,nObs,1).array() = 1/sig2*yXbeZg2.transpose().array()-1;
  // std::cout << "tG = \n" << tG << std::endl;
  // quantile param(s)
  for (int ii=0; ii<nObs; ii++) {
    for (int jj=0; jj<nQts; jj++) {
      tG.block(ii,nBet+nGam+1+jj,1,1).array() = phi_tau_smooth(yXbeZg(ii)/sqrt(sig2)-Nu(jj),tau[jj],s);
    }
  }
  G = tG.transpose();
}

#endif