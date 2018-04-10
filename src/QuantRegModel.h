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
  double alpha; 
  int nObs, nEqs, nBet, nGam;
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
  void evalG(const Ref<const VectorXd>& beta);
  void evalG(const Ref<const VectorXd>& beta, 
             const Ref<const VectorXd>& gamma);
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
  alpha = *(double*)(params);
  nObs = y.size();
  nBet = X.rows(); // X gets passed as nBet x nObs matrix
  nEqs = nBet; 
  G = MatrixXd::Zero(nEqs,nObs);
}

// setData (location-scale model) (with default ctor)
inline void QuantRegModel::setData(const Ref<const VectorXd>& _y,
                                   const Ref<const MatrixXd>& _X,
                                   const Ref<const MatrixXd>& _Z,
                                   void* params) {
  y = _y;
  X = _X;
  Z = _Z;
  alpha = *(double*)(params);
  nObs = y.size();
  nBet = X.rows(); // X gets passed as nBet x nObs matrix
  nGam = Z.rows(); // Z gets passed as nGam x nObs matrix
  nEqs = nBet+nGam+1; 
  G = MatrixXd::Zero(nEqs,nObs);
}

// revised L1 loss function for quantile regression
inline double QuantRegModel::rho_alpha(double u, double alpha) {
    return(u * (alpha - (u <= 0)));
}

// 1st derivative of rho_alpha
inline double QuantRegModel::phi_alpha(double u, double alpha) {
    return((u <= 0) - alpha);
}

// form the G matrix (location model)
inline void QuantRegModel::evalG(const Ref<const VectorXd>& beta) {
    for(int ii=0; ii<y.size(); ii++) {
        this->G.col(ii) = phi_alpha(y(ii)-X.col(ii).transpose()*beta, alpha)*X.col(ii);
    }
}

// form the G matrix (location-scale model)
inline void QuantRegModel::evalG(const Ref<const VectorXd>& beta, 
                                 const Ref<const VectorXd>& gamma) {
  eZg.array() = (-gamma.transpose()*Z).array().exp();
  yXbeZg.array() = (y.transpose()-beta.transpose()*X).array()*eZg.array();
  double phivalue; 
  for(int ii=0; ii<y.size(); ii++) {
    phivalue = phi_alpha(yXbeZg(ii), alpha);
    this->G.block(0,ii,nBet,1) = phivalue*eZg(ii)*X.col(ii);
    this->G.block(nBet,ii,nGam,1) = phivalue*yXbeZg(ii)*Z.col(ii);
  }
  G.bottomRows(1).array() = yXbeZg.array()*yXbeZg.array()-1;
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