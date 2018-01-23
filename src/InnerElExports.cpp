// port some of InnerEL to R

#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends("RcppEigen")]]
#include <RcppEigen.h>
using namespace Eigen;
#include "MeanReg.h"

// G: m x N matrix
// lambda0: m-vector of starting values
// [[Rcpp::export(".lambdaNR")]]
Rcpp::List lambdaNR(Eigen::MatrixXd G, Eigen::VectorXd lambda0,
		    int maxIter, double eps, bool verbose) {
  int nIter;
  double maxErr;
  VectorXd lambda(G.rows());
  bool not_conv;
  MeanReg IL(G.cols(), G.rows(), lambda0); // constructor
  IL.G = G; // assignment
  // IL.lambdaOld = lambda0; // NEW: removed since it's assigned by ctor
  IL.LambdaNR(nIter, maxErr, maxIter, eps);
  lambda = IL.lambdaNew; // output
  // check convergence
  not_conv = nIter >= maxIter && maxErr >= eps;
  if(verbose) {
    Rprintf("nIter = %i, maxErr = %f\n", nIter, maxErr);
  }
  return Rcpp::List::create(_["lambda"] = lambda,
			    _["convergence"] = !not_conv);
}

// Not needed
// // [[Rcpp::export(".PostSample")]]
// Eigen::MatrixXd PostSample(Eigen::MatrixXd G,
//                              int nObs, int nEqs, Eigen::VectorXd y, Eigen::MatrixXd X, 
//                              Eigen::VectorXd lambda0, 
//                              int nsamples, int nburn, Eigen::VectorXd betaInit,
//                              Eigen::VectorXd sigs, int maxIter = 100, double eps = 1e-7)  {
//   InnerEL IL(nObs, nEqs, lambda0);
//   IL.G = G; // assignment
//   Eigen::MatrixXd beta_chain = IL.PostSample(nsamples, nburn, y, X, betaInit,
//                                                sigs, maxIter, eps);
//   return(beta_chain);
// }
