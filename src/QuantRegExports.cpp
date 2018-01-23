// port some of QuantReg to R

#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends("RcppEigen")]]
#include <RcppEigen.h>
using namespace Eigen;
#include "QuantReg.h"


// [[Rcpp::export(".logELquant")]]
double logELquant(int nObs, int nEqs, Eigen::VectorXd y, Eigen::MatrixXd X,
                 Eigen::VectorXd beta, double alpha,
                 Eigen::VectorXd lambda0, int maxIter = 100, double eps = 1e-7) {
      QuantReg QR(nObs, nEqs, alpha, lambda0);
      double logELquant = QR.logEL(y,X,beta,maxIter,eps);
      return(logELquant);
}

// [[Rcpp::export(".QuantReg_post")]]
Eigen::MatrixXd QuantReg_post(int nObs, int nEqs, Eigen::VectorXd y, Eigen::MatrixXd X, 
                             double alpha, Eigen::VectorXd lambda0, 
                             int nsamples, int nburn, Eigen::VectorXd betaInit,
                             Eigen::VectorXd sigs, int maxIter = 100, double eps = 1e-7) {
      QuantReg QR(nObs, nEqs, alpha, lambda0);
      Eigen::MatrixXd beta_chain = QR.PostSample(nsamples, nburn, y, X, betaInit,
                                                   sigs, maxIter, eps);
      return(beta_chain);
}