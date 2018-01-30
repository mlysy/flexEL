// port some of InnerEL to R

#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends("RcppEigen")]]
#include <RcppEigen.h>
using namespace Eigen;
#include "InnerEL.h"
#include "MeanRegModel.h"

// Note: y, X are not actually needed here but instantiating an InnerEL object needs them 
// G: m x N matrix
// lambda0: m-vector of starting values
// [[Rcpp::export(".lambdaNR")]]
Rcpp::List lambdaNR(Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::MatrixXd G, 
                    Eigen::VectorXd lambda0, int nObs, int nEqs, 
                    int maxIter, double eps, bool verbose) {
  InnerEL<MeanRegModel> IL(y, X, nObs, nEqs, lambda0); // instantiate
  IL.G = G; // assign a given G
  // initialize variables for output here 
  int nIter;
  double maxErr;
  VectorXd lambda(nEqs);
  bool not_conv;
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