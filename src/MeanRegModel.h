#ifndef MEANREGMODEL_h
#define MEANREGMODEL_h

#include <math.h>
#include <Rmath.h>

// subclass: Mean regression
class MeanRegModel {
protected:
  VectorXd y;
  MatrixXd X;
  int nObs, nEqs;
  MatrixXd GGt;
public:
  MatrixXd G; // TODO: can make a set function instead of making it public 
  MeanRegModel(VectorXd _y, MatrixXd _X, int nObs, int nEqs); // constructor non-censor
  void evalG(VectorXd beta);
};

// constructor
inline MeanRegModel::MeanRegModel(VectorXd _y, MatrixXd _X, int nObs, int nEqs) {
  this->y = _y;
  this->X = _X; 
  this->nObs = nObs;
  this->nEqs = nEqs;
  this->G = MatrixXd::Zero(nEqs,nObs);
  this->GGt = MatrixXd::Zero(nEqs,nObs*nEqs);
}

// Chen-Van Keilegom (2009) way of constructing G, 
// but initialization of G should be of dim p x n
// double MeanRegModel::logEL(VectorXd y, MatrixXd X, VectorXd beta,
//                       int maxIter, double eps) {
//       int nIter;
//       double maxErr;
//       for (int ii=0; ii<y.size(); ii++) {
//         this->G(0,ii) = (y(ii)-X.col(ii).transpose()*beta)*X.col(ii)
//       }
//       InnerEL::LambdaNR(nIter, maxErr, maxIter, eps);
//       VectorXd logomegahat = log(1/(1-(lambdaNew.transpose()*G).array())) - 
//         log((1/(1-(lambdaNew.transpose()*G).array())).sum());
//       return(logomegahat.sum());
// }

// form the G matrix
inline void MeanRegModel::evalG(VectorXd beta) {
  for (int ii=0; ii<y.size(); ii++) {
    this->G(0,ii) = y(ii)-X.col(ii).transpose()*beta;
    this->G(1,ii) = pow((y(ii)-X.col(ii).transpose()*beta),2.0) - 1;
  }
}


#endif