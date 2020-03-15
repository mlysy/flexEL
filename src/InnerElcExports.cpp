/**
 * @file InnerElcExports.cpp
 * 
 * @brief Export InnerELC functions to R.
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include "InnerELC.h"
// #include "MeanRegModel.h"

//[[Rcpp::depends("RcppEigen")]]

using namespace Rcpp;
using namespace Eigen;

/*
// Note: for testing purpose only
// [[Rcpp::export(".EvalEpsilonsLS")]]
Eigen::VectorXd EvalEpsilonsLS(Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::MatrixXd Z,
                               Eigen::VectorXd beta, Eigen::VectorXd gamma, double sig2) {
  // InnerELC<MeanRegModel> ILC;
  int n_obs = G.cols();
  int n_eqs = G.rows();
  InnerELC<MeanRegModel> ILC(n_obs,n_eqs);
  Eigen::VectorXd deltas = VectorXd::Zero(y.size()).array()+1.0;
  // ILC.setData(y,X,Z,deltas,NULL);
  ILC.setData(y,X,Z,deltas);
  ILC.EvalEpsilons(beta,gamma,sig2);
  Eigen::VectorXd epsilons = ILC.getEpsilons();
  return(epsilons);
}
*/

// returns weights for the weighted maximum log EL
// Note: for testing purpose only
// [[Rcpp::export(".EvalWeights")]]
Eigen::VectorXd EvalWeights(Eigen::VectorXd deltas, Eigen::VectorXd omegas, 
                            Eigen::VectorXd epsilons, bool support) {
  int n_obs = epsilons.size(); // TODO: this is problematic -- which n_obs to use
  el::InnerELC ILC(n_obs,1);
  ILC.set_opts(support);
  ILC.set_deltas(deltas);
  ILC.set_omegas(omegas);
  ILC.ste_epsilons(epsilons); 
  ILC.EvalWeights(); 
  Eigen::VectorXd weights = ILC.get_weights(); 
  return(weights);
}

// [[Rcpp::export(".LambdaNRC")]]
Eigen::VectorXd LambdaNRC(Eigen::MatrixXd G, Eigen::VectorXd weights, 
                          int max_iter, double rel_tol, bool support, bool verbose) {
  int n_obs = G.cols();
  int n_eqs = G.rows();
  el::InnerELC ILC(n_obs,n_eqs);
  ILC.set_opts(max_iter, rel_tol, support);
  ILC.set_G(G); // assign a given G
  ILC.set_weights(weights);
  // ILC.setTol(max_iter, rel_tol);
  
  // initialize variables for output here 
  int n_iter;
  double max_err;
  ILC.LambdaNR(n_iter, max_err);
  VectorXd lambda = ILC.get_lambda(); // output

  // check convergence
  bool not_conv = (n_iter == max_iter) && (max_err > rel_tol);
  if(verbose) {
    Rprintf("n_iter = %i, max_err = %f\n", n_iter, max_err);
  }

  // fill in NaN if not converged
  if (not_conv) {
    for (int ii=0; ii<lambda.size(); ii++) {
      lambda(ii) = std::numeric_limits<double>::quiet_NaN();
    }
  }
  return lambda;
  
  // return the status and value
  // return Rcpp::List::create(_["lambda"] = lambda,
  //                           _["convergence"] = !not_conv);
}

// G: m x N matrix
// lambda0: m-vector of starting values
// [[Rcpp::export(".OmegaHatEM")]]
Eigen::VectorXd OmegaHatEM(Eigen::VectorXd omegas_init, 
                           Eigen::MatrixXd G, Eigen::VectorXd deltas,
                           Eigen::VectorXd epsilons, 
                           int max_iter, double rel_tol, double abs_tol, bool support, bool verbose) {
  int n_obs = G.cols();
  int n_eqs = G.rows();
  el::InnerELC ILC(n_obs,n_eqs); 
  ILC.set_opts(max_iter, rel_tol, abs_tol, support);
  ILC.set_deltas(deltas);
  ILC.set_G(G); // assign a given G
  ILC.ste_epsilons(epsilons);
  // ILC.setTol(max_iter, rel_tol, abs_tol);
  ILC.set_omegas(omegas_init); // set initial omegas from uncensored omega.hat
  ILC.EvalOmegas();
  VectorXd omegasnew = ILC.get_omegas(); // output
  return omegasnew;
  // return Eigen::VectorXd::Zero(n_obs);
}

// [[Rcpp::export(".LogELC")]]
double LogELC(Eigen::VectorXd omegas, Eigen::VectorXd epsilons, 
              Eigen::VectorXd deltas, bool support) {
  int n_obs = omegas.size();
  el::InnerELC ILC(n_obs,1);
  ILC.set_opts(support);
  ILC.set_deltas(deltas);
  ILC.ste_epsilons(epsilons); 
  ILC.set_omegas(omegas);
  double logel = ILC.logEL();
  return logel; 
}

// // [[Rcpp::export(".EvalPsos.smooth")]]
// double EvalPsosSmooth(int ii, Eigen::VectorXd omegas, 
//                       Eigen::VectorXd epsilons, double s, bool support) {
//   int n_obs = epsilons.size();
//   el::InnerELC<MeanRegModel> ILC(n_obs,1);
//   ILC.set_opts(support);
//   ILC.ste_epsilons(epsilons);
//   ILC.set_omegas(omegas);
//   return ILC.EvalPsosSmooth(ii-1,s); // ii is 1 larger in R than in C++
// }

// [[Rcpp::export(".LogELSmooth")]]
double LogELSmooth(Eigen::VectorXd omegas, 
                   Eigen::VectorXd epsilons, 
                   Eigen::VectorXd deltas, double s, bool support) {
  int n_obs = omegas.size();
  el::InnerELC ILC(n_obs,1);
  ILC.set_opts(support);
  ILC.set_omegas(omegas);
  ILC.ste_epsilons(epsilons);
  ILC.set_deltas(deltas);
  return ILC.LogELSmooth(s);
}

// [[Rcpp::export(".EvalWeightsSmooth")]]
Eigen::VectorXd EvalWeightsSmooth(Eigen::VectorXd deltas, 
                                  Eigen::VectorXd omegas, 
                                  Eigen::VectorXd epsilons, double s, bool support) {
  int n_obs = epsilons.size();
  el::InnerELC ILC(n_obs,1);
  ILC.set_opts(support);
  ILC.set_omegas(omegas);
  ILC.ste_epsilons(epsilons);
  ILC.set_deltas(deltas);
  ILC.EvalWeightsSmooth(s);
  Eigen::VectorXd weights = ILC.get_weights(); 
  return(weights);
}

// [[Rcpp::export(".OmegaHatEMSmoothh")]]
Eigen::VectorXd OmegaHatEMSmooth(Eigen::VectorXd omegas_init,
                                 Eigen::MatrixXd G, Eigen::VectorXd deltas,
                                 Eigen::VectorXd epsilons, double s, 
                                 int max_iter, double rel_tol, double abs_tol, 
                                 bool support, bool verbose) {
  int n_obs = G.cols();
  int n_eqs = G.rows();
  el::InnerELC ILC(n_obs,n_eqs); 
  ILC.set_opts(max_iter, rel_tol, abs_tol, support);
  ILC.set_deltas(deltas);
  ILC.set_G(G); // assign a given G
  ILC.ste_epsilons(epsilons); 
  // ILC.setTol(max_iter, rel_tol, abs_tol);
  ILC.set_omegas(omegas_init); // set initial omegas from uncensored omega.hat
  ILC.EvalOmegasSmooth(s);
  VectorXd omegasnew = ILC.get_omegas(); // output
  return omegasnew; 
}
