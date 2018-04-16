#ifndef QUANTREGMODEL_h
#define QUANTREGMODEL_h

// #include <math.h>
// #include <Rmath.h>

class QuantRegModel {
private:
  RowVectorXd yXb;
  RowVectorXd eZg;
  RowVectorXd yXbeZg;
  double rho_alpha(double u, double alpha); // TODO: not needed?
  double phi_alpha(double u, double alpha); 
protected:
  VectorXd y;
  MatrixXd X;
  MatrixXd Z;
  // double alpha; 
  double *alpha; // quantile levels
  int nObs, nEqs, nBet, nGam, nQts; // nQts := number of quantile levels
  MatrixXd G;
public:
  // QuantRegModel(const Ref<const VectorXd>& _y, const Ref<const MatrixXd>& _X,
  //               void* params); // constructor non-censor
  void setData(const Ref<const VectorXd>& _y, 
               const Ref<const MatrixXd>& _X,
               void* params); // set data with default ctor
  void setData(const Ref<const VectorXd>& _y, 
               const Ref<const MatrixXd>& _X,
               const Ref<const MatrixXd>& _Z,
               void* params); // set data with default ctor
  // void evalG(const Ref<const VectorXd>& beta);
  void evalG(const Ref<const MatrixXd>& Beta);
  // void evalG(const Ref<const VectorXd>& beta,
  //            const Ref<const VectorXd>& gamma);
  void evalG(const Ref<const MatrixXd>& Beta, 
             const Ref<const MatrixXd>& Gamma);
  void setG(const Ref<const MatrixXd>& _G); 
  MatrixXd getG(); // TODO: should prob move to EL since duplicate for MR and QR
};

/*
// constructor: old 
inline QuantRegModel::QuantRegModel(const Ref<const VectorXd>& _y,
                                    const Ref<const MatrixXd>& _X,
                                    void* params) {
    y = _y;
    X = _X;
    alpha = *(double*)(params);
    nObs = y.size();
    nEqs = X.rows(); // X gets passed as p x nObs matrix
    G = MatrixXd::Zero(nEqs,nObs);
}
*/

// setData (with default ctor)
inline void QuantRegModel::setData(const Ref<const VectorXd>& _y,
                                   const Ref<const MatrixXd>& _X,
                                   void* params) {
  y = _y;
  X = _X;
  // alpha = *(double*)(params);
  double *alphaArr = (double*)(params);
  nQts = alphaArr[0]; // first entry must be the number of quantile levels
  alpha = &(alphaArr[1]); // rest are the quantile levels
  
  nObs = y.size();
  nBet = X.rows(); // X gets passed as nBet x nObs matrix
  nEqs = nBet; 
  
  // G = MatrixXd::Zero(nEqs,nObs);
  G = MatrixXd::Zero(nQts*nEqs,nObs);
}

// setData (location-scale model) (with default ctor)
inline void QuantRegModel::setData(const Ref<const VectorXd>& _y,
                                   const Ref<const MatrixXd>& _X,
                                   const Ref<const MatrixXd>& _Z,
                                   void* params) {
  y = _y;
  X = _X;
  Z = _Z;
  // alpha = *(double*)(params);
  double *alphaArr = (double*)(params);
  nQts = alphaArr[0]; // first entry must be the number of quantile levels
  alpha = &(alphaArr[1]); // rest are the quantile levels
  
  nObs = y.size();
  nBet = X.rows(); // X gets passed as nBet x nObs matrix
  nGam = Z.rows(); // Z gets passed as nGam x nObs matrix
  nEqs = nBet+nGam+1; 
  
  // G = MatrixXd::Zero(nEqs,nObs);
  G = MatrixXd::Zero(nQts*nEqs,nObs);
}

// revised L1 loss function for quantile regression
inline double QuantRegModel::rho_alpha(double u, double alpha) {
    return(u * (alpha - (u <= 0)));
}

// 1st derivative of rho_alpha
inline double QuantRegModel::phi_alpha(double u, double alpha) {
    return((u <= 0) - alpha);
}

// // single quantile case
// inline void QuantRegModel::evalG(const Ref<const VectorXd>& beta) {
//   for(int ii=0; ii<y.size(); ii++) {
//       this->G.col(ii) = phi_alpha(y(ii)-X.col(ii).transpose()*beta, alpha[0])*X.col(ii);
//   }
// }

// multiple quantile case
// form the G matrix (location model)
// inline void QuantRegModel::evalG(const Ref<const VectorXd>& beta) {
inline void QuantRegModel::evalG(const Ref<const MatrixXd>& Beta) {
  // stack multiple "Gs" for each quantile level
  for(int jj=0; jj<nQts; jj++) {
    for(int ii=0; ii<y.size(); ii++) {
      this->G.block(jj*nQts,ii,nEqs,1) = phi_alpha(y(ii)-X.col(ii).transpose()*Beta.col(jj), alpha[jj])*X.col(ii);
    }
  }
  // single quantile code:
  // for(int ii=0; ii<y.size(); ii++) {
  //     this->G.col(ii) = phi_alpha(y(ii)-X.col(ii).transpose()*beta, alpha)*X.col(ii);
  // }
}

// single quantile case
// inline void QuantRegModel::evalG(const Ref<const VectorXd>& beta,
//                                  const Ref<const VectorXd>& gamma) {
//   eZg.array() = (-gamma.transpose()*Z).array().exp();
//   yXbeZg.array() = (y.transpose()-beta.transpose()*X).array()*eZg.array();
//   double phivalue;
//   for(int ii=0; ii<y.size(); ii++) {
//     phivalue = phi_alpha(yXbeZg(ii), alpha[0]);
//     this->G.block(0,ii,nBet,1) = phivalue*eZg(ii)*X.col(ii);
//     this->G.block(nBet,ii,nGam,1) = phivalue*yXbeZg(ii)*Z.col(ii);
//   }
//   G.bottomRows(1).array() = yXbeZg.array()*yXbeZg.array()-1;
// }

// multiple quantile case
// form the G matrix (location-scale model)
inline void QuantRegModel::evalG(const Ref<const MatrixXd>& Beta,
                                 const Ref<const MatrixXd>& Gamma) {
  for (int jj=0; jj<nQts; jj++) {
    eZg.array() = (-Gamma.col(jj).transpose()*Z).array().exp();
    yXbeZg.array() = (y.transpose()-Beta.col(jj).transpose()*X).array()*eZg.array();
    double phivalue;
    for(int ii=0; ii<y.size(); ii++) {
      phivalue = phi_alpha(yXbeZg(ii), alpha[jj]);
      this->G.block(jj*nQts,ii,nBet,1) = phivalue*eZg(ii)*X.col(ii);
      this->G.block(jj*nQts+nBet,ii,nGam,1) = phivalue*yXbeZg(ii)*Z.col(ii);
    }
    G.bottomRows(1).array() = yXbeZg.array()*yXbeZg.array()-1;
  }
  // single quantile code:
  // eZg.array() = (-gamma.transpose()*Z).array().exp();
  // yXbeZg.array() = (y.transpose()-beta.transpose()*X).array()*eZg.array();
  // double phivalue;
  // for(int ii=0; ii<y.size(); ii++) {
  //   phivalue = phi_alpha(yXbeZg(ii), alpha);
  //   this->G.block(0,ii,nBet,1) = phivalue*eZg(ii)*X.col(ii);
  //   this->G.block(nBet,ii,nGam,1) = phivalue*yXbeZg(ii)*Z.col(ii);
  // }
  // G.bottomRows(1).array() = yXbeZg.array()*yXbeZg.array()-1;
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