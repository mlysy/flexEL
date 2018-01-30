#ifndef MEANREGMODEL_h
#define MEANREGMODEL_h

#include <math.h>
#include <Rmath.h>

// subclass: Mean regression
class MeanRegModel {
 private:
  RowVectorXd yXb;
  MatrixXd tG;
protected:
  VectorXd y;
  MatrixXd X;
  int nObs, nEqs;
public:
  MatrixXd G; // TODO: can make a set function instead of making it public
  MeanRegModel(const Ref<const VectorXd>& _y, const Ref<const MatrixXd>& _X,
	       void* params);
  void evalG(const Ref<const VectorXd>& beta);
};

// constructor
inline MeanRegModel::MeanRegModel(const Ref<const VectorXd>& _y,
				  const Ref<const MatrixXd>& _X,
				  void* params) {
  y = _y;
  X = _X;
  nObs = y.size();
  nEqs = X.rows(); // X gets passed as p x n matrix
  G = MatrixXd::Zero(nEqs,nObs);
  tG = MatrixXd::Zero(nObs, nEqs);
  yXb = RowVectorXd::Zero(nObs);
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
inline void MeanRegModel::evalG(const Ref<const VectorXd>& beta) {
  yXb.noalias() = y.transpose() - beta.transpose() * X;
  tG = X.transpose();
  tG.array().colwise() *= yXb.transpose().array();
  G = tG.transpose();
  // G.col(0) = y - Xb;
  // G.col(1) = 
  // for (int ii=0; ii<nObs; ii++) {
  //   G(0,ii) = y(ii)-X.col(ii).transpose()*beta;
  //   G(1,ii) = pow((y(ii)-X.col(ii).transpose()*beta),2.0) - 1;
  // }
}


#endif
