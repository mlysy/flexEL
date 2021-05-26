// [[Rcpp::depends("flexEL")]]
// [[Rcpp::depends("RcppEigen")]]

#include <Rcpp.h>
#include <RcppEigen.h>
#include <flexEL/inner_el.h>
#include <flexEL/mean_reg_model.h>

// [[Rcpp::export]]
double example_logel(Eigen::VectorXd beta,
                     Eigen::MatrixXd X,
                     Eigen::VectorXd y,
                     bool verbose) {
  flexEL::MeanRegModel MR(y, X);
  int n_obs = MR.get_n_obs();
  int n_eqs = MR.get_n_eqs();
  MatrixXd G = MatrixXd::Zero(MR.get_n_eqs(), MR.get_n_obs());
  MR.EvalG(G, beta);
  // std::cout << "G = \n" << G << std::endl;
  flexEL::GenEL GEL(n_obs, n_eqs);
  GEL.set_supp_adj(true); // turn on support correction
  double logel = GEL.logel(G);
  if(verbose) {
    int n_iter;
    double max_err;
    GEL.get_diag(n_iter, max_err);
    Rprintf("n_iter = %i, max_err = %f\n", n_iter, max_err);
  }
  return(logel);
}
