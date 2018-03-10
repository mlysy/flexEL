// port some of QuantReg to R

#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends("RcppEigen")]]
#include <RcppEigen.h>
using namespace Eigen;
#include "InnerEL.h"
#include "QuantRegModel.h"


// [[Rcpp::export(".QuantReg_logEL")]]
double QuantReg_logEL(Eigen::VectorXd y, Eigen::MatrixXd X, 
                 double alpha, Eigen::VectorXd beta, 
                 int maxIter = 100, double relTol = 1e-7) {
    InnerEL<QuantRegModel> QR(y, X, &alpha); // instantiate
    int nIter;
    double maxErr;
    double logELquant = QR.logEL(nIter, maxErr, maxIter,relTol);
    // TODO: check convergence 
    return(logELquant);
}

// // [[Rcpp::export(".PostSample_QuantReg")]]
// Eigen::MatrixXd PostSample_QuantReg(Eigen::VectorXd y, Eigen::MatrixXd X, int nObs, int nEqs, 
//                              double alpha, Eigen::VectorXd lambda0, 
//                              int nsamples, int nburn, Eigen::VectorXd betaInit,
//                              Eigen::VectorXd sigs, int maxIter = 100, double eps = 1e-7) {
//       InnerEL<QuantRegModel> QR(y, X, alpha, nObs, nEqs, lambda0); // instantiate QR(nObs, nEqs, alpha, lambda0);
//       Eigen::MatrixXd beta_chain = QR.PostSample(nsamples, nburn, betaInit, sigs, maxIter, eps);
//       return(beta_chain);
// }