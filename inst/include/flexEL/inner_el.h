/**
 * @file inner_el.h
 * 
 * @brief Inner optimization for empirial likelihood problems.
 */

#ifndef INNEREL_h
#define INNEREL_h

#include <Rcpp.h>
#include <RcppEigen.h>
#include "block_outer.h" // columnwise outer product
#include "adj_G.h" // for support correction

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

/* --------------------------------------------------------------------------- */

/**
 * @brief EL namespace
 * 
 * Wrap the exported library components into a namespace called \b el to avoid potential naming conflicts with other libraries or user-defined headers.
 */
namespace flexEL {
  
  /**
   * @file       inner_el.h
   *
   * @class      InnerEL
   *
   * @brief      A template class for empirical likelihood inner optimization calculation with fully observed responses.
   */
  class InnerEL {
    
  private:
    
    // dim of G_ matrix
    int n_obs_; // number of columns
    int n_eqs_; // number of rows
    
    // constants for LogStar
    double trunc_, aa_, bb_, cc_; //constants in LogStar functions
  
    // placeholders for LambdaNR
    MatrixXd G_; // G matrix for the estimating equations
    MatrixXd G_used_;
    VectorXd lambda0_; // initial lambda in Newton-Raphson iterations
    VectorXd lambda_old_; // old lambda in Newton-Raphson iterations
    VectorXd lambda_new_; // new lambda in Newton-Raphson iterations
    VectorXd omegas_; //  empirical distribution
    MatrixXd GGt_; //  G times G transpose
    MatrixXd GGt_used_;
    VectorXd Glambda_; // G times lambda
    ArrayXd Gl11_; // values of omegas before standardization
    VectorXd Q1_; // 1st derivative of the objective function
    MatrixXd Q2_; // 2nd derivative of the objective function
    LDLT<MatrixXd> Q2ldlt_; // robust Cholesky decomposition of Q2
    VectorXd rho_; // placeholder in LambdaNR
    VectorXd rel_err_; //relative error
    
    // tolerance values for LambdaNR and support correction
    int max_iter_; // maximum number of iterations
    double rel_tol_; // relative tolerance
    bool supp_; // whether do support correction or not
    double supp_a_; // tuning parameter for support currection

    // for support modification
    int n_obs1_; // n_obs1_ = n_obs_+1: for initial space allocation
    int n_obs2_; // n_obs2_ = n_obs_+support: used n_obs_ in calculation
    
    // maximum relative error in lambda
    double MaxRelErr();
    
    // LogStar and its derivatives for the EL dual problem
    double LogStar(double x); // A support-refined log function
    double LogStar1(double x); // First derivative of LogStar
    double LogStar2(double x); // Second derivative of LogStar
    
  public:
    
    // constructors
    InnerEL();
    InnerEL(int n_obs, int n_eqs);
    
    // set options functions
    void set_opts(const int& max_iter, const double& rel_tol,
                 const bool& supp, const double& supp_a, 
                 const Ref<const VectorXd>& lambda0);
    void set_opts(const int& max_iter, const double& rel_tol,
                 const bool& supp, 
                 const Ref<const VectorXd>& lambda0);
    void set_opts(const int& max_iter, const double& rel_tol,
                 const bool& supp);
    void set_opts(const int& max_iter, const double& rel_tol,
                 const bool& supp, const double& supp_a);
    void set_opts(const bool& supp, const double& supp_a);
    void set_opts(const bool& supp);

    // core caculations
    void LambdaNR(int& n_iter, double& max_iter); // Note: rel_tol and max_iter must be set before calling
    void EvalOmegas(); // empirical distribution
    double LogEL(); // log empirical likelihood
    void LogELGrad(double& logel, MatrixXd& dldG); // calculates log empirical likelihood and gradient
    
    // set and get functions 
    void set_lambda(const Ref<const VectorXd>& lambda); // assigned to lambdaNew
    void set_omegas(const Ref<const VectorXd>& omegas); 
    void set_G(const Ref<const MatrixXd>& G);
    VectorXd get_lambda(); 
    VectorXd get_omegas();
    MatrixXd get_G(); // TODO: is it better to return a reference?..
    Ref<MatrixXd> get_ref_G();
    
  };

} // namespace flexEL

/* --------------------------------------------------------------------------- */

// private functions

// maximum relative error in lambda
inline double flexEL::InnerEL::MaxRelErr() {
  // TODO: added for numerical stability, what is a good tolerance to use ?
  if ((lambda_new_ - lambda_old_).array().abs().maxCoeff() < 1e-10) return(0);
  rel_err_ = ((lambda_new_ - lambda_old_).array() / (lambda_new_ + lambda_old_).array()).abs();
  return(rel_err_.maxCoeff());
}

// LogStar
inline double flexEL::InnerEL::LogStar(double x) {
  if(x >= trunc_) {
    return(log(x));
  } else {
    return((aa_*x + bb_)*x + cc_);
  }
}

// d LogStar(x)/dx
inline double flexEL::InnerEL::LogStar1(double x) {
  if(x >= trunc_) {
    return(1.0/x);
  } 
  else {
    // return(aa*x + bb);
    return(2.0*aa_*x + bb_); // TODO: should be (2*aa*x + bb) ?
  }
}

// d^2 LogStar(x)/dx^2
inline double flexEL::InnerEL::LogStar2(double x) {
  if(x >= trunc_) {
    return(-1.0/(x*x));
  } else {
    // return(aa);
    return(2.0*aa_); // TODO: should be 2aa ?
  }
}

/* --------------------------------------------------------------------------- */

// public functions

// constructors
/**
 * @brief Default constructor for InnerEL.
 */
inline flexEL::InnerEL::InnerEL(){}

/**
 * @brief Constructor for InnerEL with dimensions of G matrix as inputs for memory allocation.
 * 
 * @param n_obs    Number of observations.
 * @param n_eqs    Number of estimating equations.
 */
inline flexEL::InnerEL::InnerEL(int n_obs, int n_eqs) {
  // assign internal values
  n_obs_ = n_obs;
  n_eqs_ = n_eqs;
  // support correction
  supp_ = false;
  n_obs1_ = n_obs_+1;
  n_obs2_ = n_obs_+supp_;
  // initialization of constants
  trunc_ = 1.0 / n_obs_;
  aa_ = -.5 * n_obs_*n_obs_;
  bb_ = 2.0 * n_obs_;
  cc_ = -1.5 - log(n_obs_);
  // space allocation
  omegas_ = VectorXd::Zero(n_obs1_).array() + 1.0/(double)n_obs1_; // Initialize to 1/n_obs_
  G_ = MatrixXd::Zero(n_eqs_,n_obs1_); // augmented G
  GGt_ = MatrixXd::Zero(n_eqs_,n_obs1_*n_eqs_);
  lambda0_ = VectorXd::Zero(n_eqs_); // default initial value of lambda
  lambda_old_ = VectorXd::Zero(n_eqs_); // Initialize to all 0's
  lambda_new_ = VectorXd::Zero(n_eqs_);
  Q1_ = VectorXd::Zero(n_eqs_);
  Q2_ = MatrixXd::Zero(n_eqs_,n_eqs_);
  Glambda_ = VectorXd::Zero(n_obs1_);
  Gl11_ = ArrayXd::Zero(n_obs1_);
  rho_ = VectorXd::Zero(n_obs1_);
  rel_err_ = VectorXd::Zero(n_eqs_);
  Q2ldlt_.compute(MatrixXd::Identity(n_eqs_,n_eqs_));
}

// set options functions

/**
 * @brief Set tolerance values for NR, support correction option and tuning parameter, and initial value of lambda.
 * 
 * @param max_iter    Maximum number of iterations.
 * @param rel_tol     Relative tolerance.
 * @param supp    Whether to have support correction.
 * @param supp_a       Tuning parameter for support correction (referred as "a" in chen-et-al08).
 * @param lambda0    Initial value for lambda.
 */
inline void flexEL::InnerEL::set_opts(const int& max_iter, const double& rel_tol, 
                                          const bool& supp, const double& supp_a, 
                                          const Ref<const VectorXd>& lambda0) {
  max_iter_ = max_iter;
  rel_tol_ = rel_tol;
  supp_ = supp;
  supp_a_ = supp_a;
  n_obs2_ = n_obs_+supp_;
  lambda0_ = lambda0;
}

/**
 * @brief Set tolerance values for NR, support correction option, and initial value of lambda.
 * 
 * @param max_iter    Maximum number of iterations.
 * @param rel_tol     Relative tolerance.
 * @param lambda0    Initial value for lambda.
 */
inline void flexEL::InnerEL::set_opts(const int& max_iter, const double& rel_tol, 
                                          const bool& supp, 
                                          const Ref<const VectorXd>& lambda0) {
  max_iter_ = max_iter;
  rel_tol_ = rel_tol;
  supp_ = supp;
  supp_a_ = std::max(1.0,0.5*log(n_obs_));
  n_obs2_ = n_obs_+supp_;
  lambda0_ = lambda0;
}

/**
 * @brief Set tolerance values for NR, support correction option. Initial value of lambda is omitted and is default to be a vector of zeros.
 * 
 * @param max_iter    Maximum number of iterations.
 * @param rel_tol     Relative tolerance.
 * @param supp    Whether to have support correction.
 * @param supp_a       Tuning parameter for support correction (referred as "a" in chen-et-al08).
 */
inline void flexEL::InnerEL::set_opts(const int& max_iter, const double& rel_tol, 
                                          const bool& supp, const double& supp_a) {
  max_iter_ = max_iter;
  rel_tol_ = rel_tol;
  supp_ = supp;
  supp_a_ = supp_a;
  n_obs2_ = n_obs_+supp_;
}

/**
 * @brief Set tolerance values for NR, support correction option. Initial value of lambda is omitted and is default to be a vector of zeros.
 * 
 * @param max_iter    Maximum number of iterations.
 * @param rel_tol     Relative tolerance.
 * @param supp    Whether to have support correction.
 */
inline void flexEL::InnerEL::set_opts(const int& max_iter, const double& rel_tol, 
                                          const bool& supp) {
  max_iter_ = max_iter;
  rel_tol_ = rel_tol;
  supp_ = supp;
  supp_a_ = std::max(1.0,0.5*log(n_obs_));
  n_obs2_ = n_obs_+supp_;
}

/**
 * @brief Set support correction option and tuning parameter only. 
 * 
 * @param supp    Whether to have support correction.
 * @param supp_a       Tuning parameter for support correction (referred as "a" in chen-et-al08).
 */
inline void flexEL::InnerEL::set_opts(const bool& supp, const double& supp_a) {
  supp_ = supp;
  supp_a_ = supp_a;
  n_obs2_ = n_obs_+supp_;
}

/**
 * @brief Set support correction option only. 
 * 
 * @param supp    Whether to have support correction.
 */
inline void flexEL::InnerEL::set_opts(const bool& supp) {
  supp_ = supp;
  supp_a_ = std::max(1.0,0.5*log(n_obs_));
  n_obs2_ = n_obs_+supp_;
}

/**
 * @brief Find the optimal lambda by a Newton-Raphson algorithm.
 * 
 * @param[out] n_iter    Number of iterations to achieve convergence.
 * @param[out] max_iter   Maximum relative error among entires in lambda at the last step.
 */
inline void flexEL::InnerEL::LambdaNR(int& n_iter, double& max_iter) {
  
  lambda_old_ = lambda0_; // set to initial value
  lambda_new_.fill(0.0); // may not be needed here..
  
  // Note: these two cannot be preallocate untill `supp` is set
  GGt_used_ = GGt_.block(0,0,n_eqs_,n_obs2_*n_eqs_);
  G_used_ = G_.block(0,0,n_eqs_,n_obs2_);
  
  block_outer(GGt_used_,G_used_);
  
  // newton-raphson loop
  int ii, jj;
  for(ii=0; ii<max_iter_; ii++) {
    // Q1 and Q2
    Glambda_.noalias() = lambda_old_.transpose() * G_used_;
    Glambda_= 1.0 - Glambda_.array();
    Q2_.fill(0.0);
    for(jj=0; jj<n_obs2_; jj++) {
      rho_(jj) = LogStar1(Glambda_(jj));
      Q2_ += LogStar2(Glambda_(jj)) * GGt_used_.block(0,jj*n_eqs_,n_eqs_,n_eqs_);
    }
    // update lambda
    Q1_ = -G_used_ * rho_;
    Q2ldlt_.compute(Q2_);
    lambda_new_.noalias() = lambda_old_ - Q2ldlt_.solve(Q1_);
    max_iter = MaxRelErr();
    if (max_iter < rel_tol_) {
      break;
    }
    lambda_old_ = lambda_new_; // complete cycle
  }
  n_iter = ii; // output lambda and also n_iter and max_iter
  return;
}

/**
 * @brief Evaluate omegas based on G and lambdaNew.
 */
inline void flexEL::InnerEL::EvalOmegas() {
  // G and lambdaNew must have been assigned
  if (lambda_new_ != lambda_new_) { // if lambdaNew is NaN 
    for (int ii=0; ii<n_obs2_; ii++) {
      omegas_(ii) = std::numeric_limits<double>::quiet_NaN();
    }
  }
  else {
    // if (supp_) adj_G(G_,supp_a_); // calculate based on adjusted G
    MatrixXd G_used_ = G_.block(0,0,n_eqs_,n_obs2_); // TODO: new allocation right now..
    Glambda_.head(n_obs2_).noalias() = lambda_new_.transpose() * G_used_;
    Gl11_.head(n_obs2_) = 1.0/(1.0-Glambda_.head(n_obs2_).array());
    omegas_.head(n_obs2_).array() = Gl11_.head(n_obs2_).array() / Gl11_.head(n_obs2_).sum(); // in fact, Gl11.sum() should be equal to n_obs
  }
}

/**
 * @brief Calculate LogEL using omegas.
 * 
 * @return log empirical likelihood.
 */
inline double flexEL::InnerEL::LogEL() {
  // if omegas are NaN, return -Inf
  if (omegas_.head(n_obs2_) != omegas_.head(n_obs2_)) return -INFINITY;
  else return(omegas_.head(n_obs2_).array().log().sum());
}

/**
 * @brief Calculate LogEL and gradient matrix dldG
 */
inline void flexEL::InnerEL::LogELGrad(double& logel, MatrixXd& dldG) {
  if (supp_ == false) {
    if (omegas_.head(n_obs2_) != omegas_.head(n_obs2_)) {
      logel = -INFINITY;
      dldG = MatrixXd::Zero(n_eqs_,n_obs2_);
    }
    else {
      logel = omegas_.head(n_obs2_).array().log().sum();
      // std::cout << lambda_new_.transpose() << std::endl;
      // std::cout << omegas_.head(n_obs2_).transpose() << std::endl;
      dldG = lambda_new_ * omegas_.head(n_obs2_).transpose();
    }
  }
  // TODO: with support correction
}

// setters

/**
 * @brief Set the value of lambda (e.g. to be used directly to calculate omegas).
 */
inline void flexEL::InnerEL::set_lambda(const Ref<const VectorXd>& lambda) {
  // lambda_old_ = lambda;
  lambda_new_ = lambda;
}

/**
 * @brief Set the value of omegas (e.g. to be used directly to calculate log EL).
 */
inline void flexEL::InnerEL::set_omegas(const Ref<const VectorXd>& omegas) {
  omegas_.head(n_obs2_) = omegas; 
}

/**
 * @brief Set the value of G (e.g. to be used directly to calculate lambda or log EL).
 */
inline void flexEL::InnerEL::set_G(const Ref<const MatrixXd>& G) {
  G_.block(0,0,n_eqs_,n_obs_) = G;
  if (supp_) adj_G(G_,supp_a_);
}

// getters

/**
 * @brief Get the value of lambda.
 */
inline VectorXd flexEL::InnerEL::get_lambda() {
  return(lambda_new_);
}

/**
 * @brief Get the value of omegas.
 */
inline VectorXd flexEL::InnerEL::get_omegas() {
  return(omegas_.head(n_obs2_));
}

/**
 * @brief Get the value of G.
 */
inline MatrixXd flexEL::InnerEL::get_G() {
  return(G_.block(0,0,n_eqs_,n_obs2_));
}

/**
 * @brief Get the reference of original dimension G (not including the adjusted row).
 */
inline Ref<MatrixXd> flexEL::InnerEL::get_ref_G() {
  return Ref<MatrixXd>(G_.block(0,0,n_eqs_,n_obs2_-supp_));
}

#endif
