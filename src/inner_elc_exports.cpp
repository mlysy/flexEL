/**
 * @file InnerElcExports.cpp
 * 
 * @brief Export InnerELC functions to R.
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include "inner_elc.h"

//[[Rcpp::depends("RcppEigen")]]

using namespace Rcpp;
using namespace Eigen;

// returns weights for the weighted maximum log EL
// Note: for testing purpose only
// [[Rcpp::export(".EvalWeights")]]
Eigen::VectorXd EvalWeights(Eigen::VectorXd deltas, Eigen::VectorXd omegas, 
                            Eigen::VectorXd epsilons, bool support) {
  int n_obs = epsilons.size(); // TODO: this is problematic -- which n_obs to use
  flexEL::InnerELC ILC(n_obs,1);
  ILC.set_opts(support);
  ILC.set_deltas(deltas);
  ILC.set_omegas(omegas);
  ILC.set_epsilons(epsilons); 
  ILC.EvalWeights(); 
  Eigen::VectorXd weights = ILC.get_weights(); 
  return(weights);
}

// [[Rcpp::export(".LambdaNRC")]]
Eigen::VectorXd LambdaNRC(Eigen::MatrixXd G, Eigen::VectorXd weights, 
                          int max_iter, double rel_tol, bool support, bool verbose) {
  int n_obs = G.cols();
  int n_eqs = G.rows();
  flexEL::InnerELC ILC(n_obs,n_eqs);
  ILC.set_opts(max_iter, rel_tol, support);
  ILC.set_G(G); // assign a given G
  ILC.set_weights(weights);

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
  flexEL::InnerELC ILC(n_obs,n_eqs); 
  ILC.set_opts(max_iter, rel_tol, abs_tol, support);
  ILC.set_deltas(deltas);
  ILC.set_G(G); // assign a given G
  ILC.set_epsilons(epsilons);
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
  flexEL::InnerELC ILC(n_obs,1);
  ILC.set_opts(support);
  ILC.set_deltas(deltas);
  ILC.set_epsilons(epsilons); 
  ILC.set_omegas(omegas);
  double logel = ILC.LogEL();
  return logel; 
}

// [[Rcpp::export(".LogELSmooth")]]
double LogELSmooth(Eigen::VectorXd omegas, 
                   Eigen::VectorXd epsilons, 
                   Eigen::VectorXd deltas, double s, bool support) {
  int n_obs = omegas.size();
  flexEL::InnerELC ILC(n_obs,1);
  ILC.set_opts(support);
  ILC.set_omegas(omegas);
  ILC.set_epsilons(epsilons);
  ILC.set_deltas(deltas);
  return ILC.LogELSmooth(s);
}

// [[Rcpp::export(".EvalWeightsSmooth")]]
Eigen::VectorXd EvalWeightsSmooth(Eigen::VectorXd deltas, 
                                  Eigen::VectorXd omegas, 
                                  Eigen::VectorXd epsilons, double s, bool support) {
  int n_obs = epsilons.size();
  flexEL::InnerELC ILC(n_obs,1);
  ILC.set_opts(support);
  ILC.set_omegas(omegas);
  ILC.set_epsilons(epsilons);
  ILC.set_deltas(deltas);
  ILC.EvalWeightsSmooth(s);
  Eigen::VectorXd weights = ILC.get_weights(); 
  return(weights);
}

// [[Rcpp::export(".OmegaHatEMSmooth")]]
Eigen::VectorXd OmegaHatEMSmooth(Eigen::VectorXd omegas_init,
                                 Eigen::MatrixXd G, Eigen::VectorXd deltas,
                                 Eigen::VectorXd epsilons, double s, 
                                 int max_iter, double rel_tol, double abs_tol, 
                                 bool support, bool verbose) {
  int n_obs = G.cols();
  int n_eqs = G.rows();
  flexEL::InnerELC ILC(n_obs,n_eqs); 
  ILC.set_opts(max_iter, rel_tol, abs_tol, support);
  ILC.set_deltas(deltas);
  ILC.set_G(G); // assign a given G
  ILC.set_epsilons(epsilons); 
  // ILC.setTol(max_iter, rel_tol, abs_tol);
  ILC.set_omegas(omegas_init); // set initial omegas from uncensored omega.hat
  ILC.EvalOmegasSmooth(s);
  VectorXd omegasnew = ILC.get_omegas(); // output
  return omegasnew; 
}
