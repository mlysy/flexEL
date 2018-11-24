#ifndef QUANTREGMODEL_h
#define QUANTREGMODEL_h

// #include <math.h>
// #include <Rmath.h>
#include "IndSmooth.h"

class QuantRegModel {
private:
  RowVectorXd yXb;
  RowVectorXd eZg;
  RowVectorXd yXbeZg;
  RowVectorXd yXbeZg2;
  // double rho_alpha(double u, double alpha); // TODO: not needed?
  double phi_alpha(double u, double alpha); 
  MatrixXd tG;
protected:
  VectorXd y;
  MatrixXd X;
  MatrixXd Z;
  double *alpha; // quantile levels
  int nObs, nEqs, nBet, nGam, nQts; // nQts := number of quantile levels
  MatrixXd G;
public:
  QuantRegModel(); // default ctor -- it shouldn't have one actually
  QuantRegModel(int _nObs, int _nEqs);
  void setData(const Ref<const VectorXd>& _y, 
               const Ref<const MatrixXd>& _X,
               void* params); // set data with default ctor
  void setData(const Ref<const VectorXd>& _y, 
               const Ref<const MatrixXd>& _X,
               const Ref<const MatrixXd>& _Z,
               void* params); // set data with default ctor
  void evalG(const Ref<const MatrixXd>& Beta);
  void evalG(const Ref<const VectorXd>& beta, 
             const Ref<const VectorXd>& gamma,
             const double& sig2, // sig2 should be a scalar
             const Ref<const VectorXd>& Nu);
  void evalGSmooth(const Ref<const VectorXd>& beta,
                   const Ref<const VectorXd>& gamma,
                   const double& sig2, // sig2 should be a scalar
                   const Ref<const VectorXd>& Nu,
                   const double s);
  void setG(const Ref<const MatrixXd>& _G); 
  MatrixXd getG(); // TODO: should prob move to EL since duplicate for MR and QR
  // smooth functions
  // double rho_alpha_smooth(double u, double alpha);
  double phi_alpha_smooth(double u, double alpha, double s); 
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
  double *alphaArr = (double*)(params);
  nQts = alphaArr[0]; // first entry must be the number of quantile levels
  alpha = &(alphaArr[1]); // rest are the quantile levels
  
  // std::cout << "nQts = " << nQts << std::endl;
  // std::cout << "alpha = ";
  // for (int ii=0; ii<nQts; ii++) {
  //   std::cout << alpha[ii] << ' ';
  // }
  // std::cout << std::endl;
  
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
  double *alphaArr = (double*)(params);
  nQts = alphaArr[0]; // first entry must be the number of quantile levels
  alpha = &(alphaArr[1]); // rest are the quantile levels
  
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

// revised L1 loss function for quantile regression
// inline double QuantRegModel::rho_alpha(double u, double alpha) {
//     return(u * (alpha - (u <= 0)));
// }

// 1st derivative of rho_alpha
inline double QuantRegModel::phi_alpha(double u, double alpha) {
    return((u <= 0) - alpha);
}

// 1st derivative of rho_alpha_smooth
inline double QuantRegModel::phi_alpha_smooth(double u, double alpha, double s) {
  return(alpha-ind_smooth(u,s)-u*ind1_smooth(u,s));
}

/*
// single quantile case (location model)
inline void QuantRegModel::evalG(const Ref<const VectorXd>& beta) {
  for(int ii=0; ii<y.size(); ii++) {
      this->G.col(ii) = phi_alpha(y(ii)-X.col(ii).transpose()*beta, alpha[0])*X.col(ii);
  }
}
*/

// multiple quantile case (location model)
inline void QuantRegModel::evalG(const Ref<const MatrixXd>& Beta) {
  // stack multiple "Gs" for each quantile level
  for(int jj=0; jj<nQts; jj++) {
    for(int ii=0; ii<y.size(); ii++) {
      this->G.block(jj*nBet,ii,nBet,1) = phi_alpha(y(ii)-X.col(ii).transpose()*Beta.col(jj), alpha[jj])*X.col(ii);
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
      tG.block(ii,nBet+nGam+1+jj,1,1).array() = phi_alpha(yXbeZg(ii)/sqrt(sig2)-Nu(jj), alpha[jj]);
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
      tG.block(ii,nBet+nGam+1+jj,1,1).array() = phi_alpha_smooth(yXbeZg(ii)/sqrt(sig2)-Nu(jj),alpha[jj],s);
    }
  }
  G = tG.transpose();
}


// set function for G matrix
inline void QuantRegModel::setG(const Ref<const MatrixXd>& _G) {
    G = _G; 
}

// get function for G matrix
inline MatrixXd QuantRegModel::getG() {
    return(G);
}

#endif