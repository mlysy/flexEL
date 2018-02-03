// port some of MeanReg functions to R

#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends("RcppEigen")]]
#include <RcppEigen.h>
using namespace Eigen;
#include "InnerEL.h"
#include "MeanRegModel.h"


// [[Rcpp::export(".MeanReg_logEL")]]
double MeanReg_logEL(Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::VectorXd beta,
                     int maxIter = 100, double eps = 1e-7) {
    InnerEL<MeanRegModel> MR(y, X, NULL); // instantiate
    double logELmean = MR.logEL(beta,maxIter,eps);
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
