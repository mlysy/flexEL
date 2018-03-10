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
Rcpp::List lambdaNR(Eigen::MatrixXd G, 
                    int maxIter, double relTol, bool verbose) {
  int nObs = G.cols();
  int nEqs = G.rows();
  VectorXd y = VectorXd::Zero(nObs);
  MatrixXd X = MatrixXd::Zero(nEqs, nObs);
  InnerEL<MeanRegModel> IL(y, X, NULL); // instantiate
  IL.G = G; // assign a given G
  // initialize variables for output here 
  int nIter;
  double maxErr;
  VectorXd lambda(nEqs);
  bool not_conv;
  IL.LambdaNR(nIter, maxErr, maxIter, relTol);
  lambda = IL.lambdaNew; // output
  // check convergence
  not_conv = (nIter == maxIter) && (maxErr > relTol);
  if(verbose) {
    Rprintf("nIter = %i, maxErr = %f\n", nIter, maxErr);
  }
  return Rcpp::List::create(_["lambda"] = lambda,
                            _["convergence"] = !not_conv);
}

// Eigen::VectorXd y, Eigen::MatrixXd X not needed here since there is no ordering 
// as in the censored case 
// [[Rcpp::export(".omega.hat")]]
Rcpp::List evalOmegas(Eigen::MatrixXd G, 
                          int maxIter, double relTol, bool verbose) {
    int nObs = G.cols();
    int nEqs = G.rows();
    VectorXd y = VectorXd::Zero(nObs);
    MatrixXd X = MatrixXd::Zero(nEqs, nObs);
    InnerEL<MeanRegModel> IL(y, X, NULL); // instantiate
    IL.G = G; // assign G before optmization 
    int nIter;
    double maxErr;
    bool not_conv;
    IL.evalOmegas(nIter, maxErr, maxIter, relTol); // calculate omegas
    VectorXd omegasnew = IL.getOmegas(); // get omegas
    // check convergence
    not_conv = (nIter == maxIter) && (maxErr > relTol);
    if(verbose) {
        Rprintf("nIter = %i, maxErr = %f\n", nIter, maxErr);
    }
    return Rcpp::List::create(_["omegas"] = omegasnew,
                              _["convergence"] = !not_conv);
}
