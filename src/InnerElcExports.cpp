// port some of InnerELC to R

#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends("RcppEigen")]]
#include <RcppEigen.h>
using namespace Eigen;
#include "InnerELC.h"
#include "MeanRegModel.h"

// Note: y, X are not actually needed here but instantiating an InnerEL object needs them 
// G: m x N matrix
// lambda0: m-vector of starting values
// [[Rcpp::export(".lambdaNRC")]]
Rcpp::List lambdaNRC(Eigen::MatrixXd G, Eigen::VectorXd weights, 
                     int maxIter, double relTol, bool verbose) {
    // pseudo input since they are not used in calculation of lambda
    int nObs = G.cols();
    int nEqs = G.rows();
    VectorXd y = VectorXd::Zero(nObs);
    MatrixXd X = MatrixXd::Zero(nEqs, nObs);
    VectorXd delta = VectorXd::Zero(nObs);
    InnerELC<MeanRegModel> IL(y, X, delta, NULL); // instantiate
    IL.G = G; // assign a given G
    IL.weights = weights; 
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
