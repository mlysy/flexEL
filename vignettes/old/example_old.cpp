// // [[Rcpp::depends("flexEL")]]
// 
// #include <Rcpp.h>
// #include <RcppEigen.h>
// #include <flexEL/inner_el.h>
// #include <flexEL/mean_reg_model.h> // using mean regression as an example, this could be your model C++ header file
// 
// // [[Rcpp::export]]
// double example_logEL(int N = 20) {
// 
//   Rcpp::NumericVector bb = Rcpp::NumericVector::create(1.0, 2.0);
//   Eigen::VectorXd bet = Rcpp::as<Eigen::Map<Eigen::VectorXd> >(bb);
//   
//   Rcpp::NumericVector ee = Rcpp::rnorm(N);
//   Eigen::VectorXd eps = Rcpp::as<Eigen::Map<Eigen::VectorXd> >(ee);
//   
//   Eigen::MatrixXd X = Eigen::MatrixXd::Random(2,N);
//   
//   Eigen::VectorXd y = bet.transpose() * X + eps;
//   
//   // instantiate the model and EL objects
//   flexEL::MeanRegModel MR(y, X); // instantiate a mean regression object
//   flexEL::InnerEL EL(MR.get_n_obs(), MR.get_n_eqs()); // instantiate an EL object with dimention of G matrix
//   bool support = true;
//   EL.set_opts(support); // set support correction
//   
//   // evaluate EL for a particular beta
//   Eigen::VectorXd beta = Eigen::VectorXd::Random(2);
//   MR.EvalG(EL.get_ref_G(), beta);
//   int n_iter;
//   double max_err;
//   EL.LambdaNR(n_iter, max_err);
//   EL.EvalOmegas();
//   return(EL.LogEL());
// }