#ifndef QUANTREGMODEL_h
#define QUANTREGMODEL_h

#include <math.h>
#include <Rmath.h>

// subclass: Mean regression
class QuantRegModel {
protected:
  VectorXd y;
  MatrixXd X;
  double alpha; 
  int nObs, nEqs;
  MatrixXd G;
  MatrixXd GGt; 
  double rho_alpha(double u, double alpha); // TODO: not needed?
  double phi_alpha(double u, double alpha); 
public:
  QuantRegModel(VectorXd _y, MatrixXd _X, double alpha, int nObs, int nEqs); // constructor non-censor
  void evalG(VectorXd beta);
};

// constructor
inline QuantRegModel::QuantRegModel(VectorXd _y, MatrixXd _X, double alpha, int nObs, int nEqs) {
  this->y = _y;
  this->X = _X; 
  this->alpha = alpha; 
  this->nObs = nObs;
  this->nEqs = nEqs;
  this->G = MatrixXd::Zero(nEqs,nObs);
  this->GGt = MatrixXd::Zero(nEqs,nObs*nEqs);
}

// revised L1 loss function for quantile regression
inline double QuantRegModel::rho_alpha(double u, double alpha) {
  return(u * (alpha - (u <= 0)));
}

// 1st derivative of rho_alpha
inline double QuantRegModel::phi_alpha(double u, double alpha) {
  return((u <= 0) - alpha);
}

// form the G matrix
inline void QuantRegModel::evalG(VectorXd beta) {
  for(int ii=0; ii<y.size(); ii++) {
    this->G.col(ii) = phi_alpha(y(ii)-X.col(ii).transpose()*beta, alpha)*X.col(ii);
  }
}

#endif