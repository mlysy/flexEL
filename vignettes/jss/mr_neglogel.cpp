// [[Rcpp::depends("flexEL")]]
// [[Rcpp::depends("RcppEigen")]]

#include <Rcpp.h>
#include <RcppEigen.h>
#include <flexEL/inner_el.h>
#include <flexEL/mean_reg_model.h>

// [[Rcpp::export]]
double mr_neglogel(Eigen::VectorXd beta,
                   Eigen::VectorXd y,
                   Eigen::MatrixXd X,
                   bool verbose) {
  flexEL::MeanRegModel MR(y, X); // initiate MeanRegModel object with data
  int n_obs = MR.get_n_obs();
  int n_eqs = MR.get_n_eqs();
  MatrixXd G = MatrixXd::Zero(n_eqs, n_obs);
  MR.EvalG(G, beta); // obtain G matrix given certain parameter values
  flexEL::GenEL GEL(n_obs, n_eqs); // initiate GenEL object given dimensions
  GEL.set_supp_adj(true); // turn on support correction
  double logel = GEL.logel(G); // obtain log EL
  if(verbose) {
    int n_iter;
    double max_err;
    GEL.get_diag(n_iter, max_err); // get diagnostics
    Rprintf("n_iter = %i, max_err = %f\n", n_iter, max_err);
  }
  return(-logel);
}
