#ifndef MEANREG_h
#define MEANREG_h

#include "InnerEL.h"
#include <math.h>
#include <Rmath.h>

// subclass: Mean regression
class MeanReg: public InnerEL {
 public:
  MeanReg(int, int, VectorXd); // constructor
  // NEW: The following two moved to InnerEL class
  void Gfun(VectorXd y, MatrixXd X, VectorXd beta); // NEW: form G
  double logEL(VectorXd y, MatrixXd X, VectorXd beta, int maxIter, double eps);
  MatrixXd MeanReg_post(int nsamples, int nburn, VectorXd y, MatrixXd X, 
			VectorXd betaInit, VectorXd sigs, int maxIter, double eps);
};

// constructor
// inline MeanReg::MeanReg(int N, int m, VectorXd lambda0) : InnerEL(N,m) {
//       this->lambdaOld = lambda0;
// }
inline MeanReg::MeanReg(int N, int m, VectorXd lambda0) : InnerEL(N,m,lambda0) {}

// Chen-Van Keilegom (2009) way of constructing G, 
// but initialization of G should be of dim p x n
// double MeanReg::logEL(VectorXd y, MatrixXd X, VectorXd beta,
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
inline void MeanReg::Gfun(VectorXd y, MatrixXd X, VectorXd beta) {
  for (int ii=0; ii<y.size(); ii++) {
    this->G(0,ii) = y(ii)-X.col(ii).transpose()*beta;
    // TODO: DO not use pow, use temp * temp!!!!
    this->G(1,ii) = pow((y(ii)-X.col(ii).transpose()*beta),2.0) - 1;
  }
}

// log empirical likelihood for linear regression Y = X'beta + eps
// beta is ncol(X) x 1 = p x 1 vector
inline double MeanReg::logEL(VectorXd y, MatrixXd X, VectorXd beta,
			     int maxIter, double eps) {
  int nIter;
  double maxErr;
  // (moved to Gfun above)
  // for (int ii=0; ii<y.size(); ii++) {
  //       this->G(0,ii) = y(ii)-X.col(ii).transpose()*beta;
  //       this->G(1,ii) = pow((y(ii)-X.col(ii).transpose()*beta),2.0) - 1;
  // }
  Gfun(y, X, beta);
  InnerEL::LambdaNR(nIter, maxErr, maxIter, eps);
  VectorXd logomegahat = log(1/(1-(lambdaNew.transpose()*G).array())) - 
    log((1/(1-(lambdaNew.transpose()*G).array())).sum());
  return(logomegahat.sum());
}

// MatrixXd MeanReg::MeanReg_post(int nsamples, int nburn,
//                               VectorXd y, MatrixXd X, VectorXd betaInit,
//                               VectorXd sigs, int maxIter, double eps) {
//       VectorXd betaOld = betaInit;
//       VectorXd betaNew = betaOld;
//       VectorXd betaProp = betaOld;
//       int betalen = betaInit.size();
//       MatrixXd beta_chain(betaInit.size(),nsamples);
//       double logELOld = logEL(y, X, betaOld, maxIter, eps); // NEW: chache old
// 
//       for (int ii=-nburn; ii<nsamples; ii++) {
//             for (int jj=0; jj<betalen; jj++) {
//                   betaProp = betaOld;
//                   betaProp(jj) = betaOld(jj) + sigs(jj)*R::norm_rand();
//                   // check if proposed beta satisfies the constraint
//                   bool satisfy = false;
//                   int nIter = 0;
//                   double maxErr;
//                   // had a BUG here?! Didn't change G!!!
//                   Gfun(y,X,betaProp); // NEW: change G with betaProp
//                   InnerEL::LambdaNR(nIter, maxErr, maxIter, eps);
//                   if (nIter < maxIter) satisfy = true;
//                   // if does not satisfy, keep the old beta
//                   if (satisfy == false) break;
//                   // if does satisfy, flip a coin
//                   double u = R::unif_rand();
//                   // use the lambda calculate just now to get the logEL for Prop
//                   // to avoid an extra call of lambdaNR
//                   VectorXd logomegahat = log(1/(1-(lambdaNew.transpose()*G).array())) - 
//                     log((1/(1-(lambdaNew.transpose()*G).array())).sum());
//                   double logELProp = logomegahat.sum();
//                   double ratio = exp(logELProp-logELOld);
//                   // double ratio = exp(logEL(y, X, betaProp, maxIter, eps) -
//                   //                    logEL(y, X, betaOld, maxIter, eps));
//                   double a = std::min(1.0,ratio);
//                   if (u < a) { // accepted
//                         betaNew = betaProp;
//                         betaOld = betaNew;
//                         logELOld = logELProp; // NEW: store the new one
//                   }
//             }
//             if (ii >= 0) {
//                   beta_chain.col(ii) = betaNew;
//             }
//       }
//       return(beta_chain);
// }

#endif
