// port some of InnerELC to R

#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends("RcppEigen")]]
#include <RcppEigen.h>
using namespace Eigen;
#include "InnerELC.h"
#include "MeanRegModel.h"


// TODO: getWeights does not need y X, only depends on deltas, epsilons and omegas, 
//       remove y X and use pesudo input to init
// returns weights for the weighted NR algorithm 
// [[Rcpp::export(".evalWeights")]]
Eigen::VectorXd evalWeights(Eigen::VectorXd deltas, Eigen::VectorXd omegas, 
                           Eigen::VectorXd epsilons) {
    // pseudo input since they are not used in calculation of lambda
    int nObs = omegas.size();
    VectorXd y = VectorXd::Zero(nObs);
    MatrixXd X = MatrixXd::Zero(1,nObs);
    InnerELC<MeanRegModel> ILC(y, X, deltas, NULL); // instantiate
    ILC.omegas = omegas;
    ILC.setEpsilons(epsilons); 
    ILC.evalWeights(); 
    Eigen::VectorXd weights = ILC.weights; 
    return(weights);
}
// Eigen::VectorXd evalWeights(Eigen::VectorXd y, Eigen::MatrixXd X, 
//                            Eigen::VectorXd deltas, Eigen::VectorXd omegas, 
//                            Eigen::VectorXd epsilons) {
//     InnerELC<MeanRegModel> ILC(y, X, deltas, NULL); // instantiate
//     ILC.omegas = omegas;
//     ILC.setEpsilons(epsilons); 
//     ILC.evalWeights(); 
//     Eigen::VectorXd weights = ILC.weights; 
//     return(weights);
// }

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
    InnerELC<MeanRegModel> ILC(y, X, delta, NULL); // instantiate
    ILC.G = G; // assign a given G
    ILC.weights = weights; 
    // initialize variables for output here 
    int nIter;
    double maxErr;
    VectorXd lambda(nEqs);
    bool not_conv;
    ILC.LambdaNR(nIter, maxErr, maxIter, relTol);
    lambda = ILC.lambdaNew; // output
    // check convergence
    not_conv = (nIter == maxIter) && (maxErr > relTol);
    if(verbose) {
        Rprintf("nIter = %i, maxErr = %f\n", nIter, maxErr);
    }
    return Rcpp::List::create(_["lambda"] = lambda,
                              _["convergence"] = !not_conv);
}

// G: m x N matrix
// lambda0: m-vector of starting values
// [[Rcpp::export(".omega.hat.EM")]]
Rcpp::List evalOmegasEM(Eigen::MatrixXd G, Eigen::VectorXd deltas,
                Eigen::VectorXd epsilons, 
                int maxIter, double relTol, bool verbose) {
    // pseudo input since they are not used in calculation of lambda
    int nObs = G.cols();
    int nEqs = G.rows();
    VectorXd y = VectorXd::Zero(nObs);
    MatrixXd X = MatrixXd::Zero(nEqs, nObs);
    InnerELC<MeanRegModel> ILC(y, X, deltas, NULL); // instantiate
    ILC.G = G; // assign a given G
    ILC.deltas = deltas; 
    ILC.setEpsilons(epsilons); 
    // initialize variables for output here 
    int nIter;
    double maxErr;
    bool not_conv;
    ILC.evalOmegas(nIter, maxErr, maxIter, relTol);
    VectorXd omegasnew = ILC.getOmegas(); // output
    // check convergence
    not_conv = (nIter == maxIter) && (maxErr > relTol);
    if(verbose) {
        Rprintf("nIter = %i, maxErr = %f\n", nIter, maxErr);
    }
    return Rcpp::List::create(_["omegas"] = omegasnew,
                              _["convergence"] = !not_conv);
}
