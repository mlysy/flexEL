// port some of QuantReg to R

#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends("RcppEigen")]]
#include <RcppEigen.h>
using namespace Eigen;
#include "InnerEL.h"
#include "QuantRegModel.h"

// [[Rcpp::export(".QuantReg_evalG")]]
Eigen::MatrixXd QuantReg_evalG(Eigen::VectorXd y, Eigen::MatrixXd X,
                              double alpha, Eigen::VectorXd beta) {
    InnerEL<QuantRegModel> QR(y, X, &alpha); // instantiate
    QR.evalG(beta);
    Eigen::MatrixXd G = QR.getG(); // G is nEqs x nObs
    return(G); 
}

// [[Rcpp::export(".QuantReg_logEL")]]
double QuantReg_logEL(Eigen::VectorXd y, Eigen::MatrixXd X, 
                 double alpha, Eigen::VectorXd beta, 
                 int maxIter = 100, double relTol = 1e-7) {
    InnerEL<QuantRegModel> QR(y, X, &alpha); // instantiate
    QR.evalG(beta); 
    int nIter;
    double maxErr;
    double logELquant = QR.logEL(nIter, maxErr, maxIter, relTol);
    // TODO: check convergence 
    return(logELquant);
}

/*
// [[Rcpp::export(".QuantReg_post")]]
Eigen::MatrixXd QuantReg_post(Eigen::VectorXd y, Eigen::MatrixXd X, 
                             double alpha, int nsamples, int nburn, 
                             Eigen::VectorXd betaInit, Eigen::VectorXd sigs, 
                             int maxIter = 100, double relTol = 1e-7) {
      InnerEL<QuantRegModel> QR(y, X, &alpha); // instantiate QR(nObs, nEqs, alpha, lambda0);
      Eigen::MatrixXd beta_chain = QR.PostSample(nsamples, nburn, betaInit, sigs, maxIter, relTol);
      return(beta_chain);
}
*/