// port some of InnerEL to R

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
Rcpp::List lambdaNRC(Eigen::MatrixXd G, Eigen::VectorXd qs, 
                     int maxIter, double eps, bool verbose) {
    int nObs = G.cols();
    int nEqs = G.rows();
    VectorXd y = VectorXd::Zero(nObs);
    MatrixXd X = MatrixXd::Zero(nEqs, nObs);
    InnerELC<MeanRegModel> IL(y, X, NULL); // instantiate
    IL.G = G; // assign a given G
    // initialize variables for output here 
    int nIter;
    double maxErr;
    VectorXd lambda(nEqs);
    bool not_conv;
    IL.LambdaNR(qs, nIter, maxErr, maxIter, eps);
    lambda = IL.lambdaNew; // output
    // check convergence
    not_conv = (nIter == maxIter) && (maxErr > eps);
    if(verbose) {
        Rprintf("nIter = %i, maxErr = %f\n", nIter, maxErr);
    }
    return Rcpp::List::create(_["lambda"] = lambda,
                              _["convergence"] = !not_conv);
}
