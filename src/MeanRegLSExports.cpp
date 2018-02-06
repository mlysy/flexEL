// port some of MeanReg functions to R

#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends("RcppEigen")]]
#include <RcppEigen.h>
using namespace Eigen;
#include "InnerEL.h"
#include "MeanRegLSModel.h"

// [[Rcpp::export(".MeanRegLS_G")]]
Eigen::MatrixXd MeanRegLS_G(Eigen::VectorXd y,
			    Eigen::MatrixXd X,
			    Eigen::MatrixXd Z,
			    Eigen::VectorXd theta) {
  InnerEL<MeanRegLSModel> MR(y, X, Z, NULL); // instantiate
  MR.evalG(theta);
  MatrixXd G = MR.G;
  return(G);
}


// [[Rcpp::export(".MeanRegLS_logEL")]]
double MeanRegLS_logEL(Eigen::VectorXd y,
		       Eigen::MatrixXd X,
		       Eigen::MatrixXd Z,
		       Eigen::VectorXd theta,
		       int maxIter = 100, double eps = 1e-7) {
  InnerEL<MeanRegLSModel> MR(y, X, Z, NULL); // instantiate
  double logELmean = MR.logEL(theta,maxIter,eps);
  return(logELmean);
}

/*
// [[Rcpp::export(".MeanReg_PostSample")]]
Eigen::MatrixXd MeanReg_PostSample(Eigen::VectorXd y, Eigen::MatrixXd X, 
                                   int nObs, int nEqs, Eigen::VectorXd lambda0, 
int nsamples, int nburn, Eigen::VectorXd betaInit,
Eigen::VectorXd sigs, int maxIter = 100, double eps = 1e-7)  {
InnerEL<MeanRegModel> MR(y, X, nObs, nEqs, lambda0); // instantiate
Eigen::MatrixXd beta_chain = MR.PostSample(nsamples, nburn, betaInit, sigs, maxIter, eps);
return(beta_chain);
}
*/
