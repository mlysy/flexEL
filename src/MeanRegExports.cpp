// port some of MeanReg to R

#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends("RcppEigen")]]
#include <RcppEigen.h>
using namespace Eigen;
#include "MeanReg.h"


// [[Rcpp::export(".logELmean")]]
double logELmean(int nObs, int nEqs, Eigen::VectorXd y, Eigen::MatrixXd X, 
                 Eigen::VectorXd beta, Eigen::VectorXd lambda0, 
                 int maxIter = 100, double eps = 1e-7) {
      MeanReg MR(nObs, nEqs, lambda0);
      double logELmean = MR.logEL(y,X,beta,maxIter,eps);
      return(logELmean);
}

// [[Rcpp::export(".MeanReg_post")]]
Eigen::MatrixXd MeanReg_post(int nObs, int nEqs, Eigen::VectorXd y, Eigen::MatrixXd X, 
                             Eigen::VectorXd lambda0, 
                             int nsamples, int nburn, Eigen::VectorXd betaInit,
                             Eigen::VectorXd sigs, int maxIter = 100, double eps = 1e-7)  {
      MeanReg MR(nObs, nEqs, lambda0);
      Eigen::MatrixXd beta_chain = MR.PostSample(nsamples, nburn, y, X, betaInit,
                                                   sigs, maxIter, eps);
      return(beta_chain);
}