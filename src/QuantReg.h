#ifndef QUANTREG_h
#define QUANTREG_h

#include <math.h> // for exp, abs
#include <random> // for random normal
#include <Rmath.h> // for random number generation
#include "InnerEL.h"

// subclass: Quantile regression
class QuantReg: public InnerEL {
private:
      double alpha; // quantile
      double rho_alpha(double u, double alpha); // TODO: not needed?
      double phi_alpha(double u, double alpha); 
public:
      // VectorXd betaOld; // params for linear regression before update
      // VectorXd betaNew; // params for linear regression after update
      // VectorXd acc_rate; // vector of acceptance rate for each beta
      // MatrixXd beta_chain; // MCMC chain
      QuantReg(int, int, double alpha, VectorXd lambda0); // constructor
      // NEW: The following two moved to InnerEL class
      void Gfun(VectorXd y, MatrixXd X, VectorXd beta); // NEW: form G
      double logEL(VectorXd y, MatrixXd X, VectorXd beta, int maxIter, double eps);
      MatrixXd QuantReg_post(int nsamples, int nburn, VectorXd y, MatrixXd X, 
                       VectorXd betaInit, VectorXd sigs, 
                       int maxIter, double eps);
};

// revised L1 loss function for quantile regression
double QuantReg::rho_alpha(double u, double alpha) {
      return(u * (alpha - (u <= 0)));
}

// 1st derivative of rho_alpha
double QuantReg::phi_alpha(double u, double alpha) {
      return((u <= 0) - alpha);
}

// constructor
// inline QuantReg::QuantReg(int N, int m, double alpha, VectorXd lambda0): InnerEL(N,m) {
//       this->alpha = alpha;
//       this->lambdaOld = lambda0;
// }
inline QuantReg::QuantReg(int N, int m, double alpha, VectorXd lambda0): InnerEL(N,m,lambda0) {
  this->alpha = alpha;
}

// form the G matrix
void QuantReg::Gfun(VectorXd y, MatrixXd X, VectorXd beta) {
      for(int ii=0; ii<y.size(); ii++) {
        this->G.col(ii) = phi_alpha(y(ii)-X.col(ii).transpose()*beta, alpha)
        *X.col(ii);
      }
}

// log empirical likelihood
double QuantReg::logEL(VectorXd y, MatrixXd X, VectorXd beta, 
                       int maxIter, double eps) {
      int nIter;
      double maxErr;
      // (moved to Gfun above)
      // for(int ii=0; ii<y.size(); ii++) {
      //       this->G.col(ii) = phi_alpha(y(ii)-X.col(ii).transpose()*beta, alpha)
      //                                                             *X.col(ii);
      // }
      Gfun(y, X, beta);
      VectorXd lambdahat(getnDims());
      InnerEL::LambdaNR(nIter, maxErr, maxIter, eps);
      VectorXd logomegahat = log(1/(1-(lambdaNew.transpose() * G).array())) -
            log((1/(1-(lambdaNew.transpose() * G).array())).sum());
      return(logomegahat.sum());
}

// posterior sampling
// MatrixXd QuantReg::QuantReg_post(int nsamples, int nburn, 
//                        VectorXd y, MatrixXd X, VectorXd betaInit,
//                        VectorXd sigs, int maxIter, double eps) {
//       VectorXd betaOld = betaInit;
//       VectorXd betaNew = betaOld;
//       VectorXd betaProp = betaOld;
//       int betalen = betaInit.size();
//       MatrixXd beta_chain(betaInit.size(), nsamples);
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
//                   if (u < a) {
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