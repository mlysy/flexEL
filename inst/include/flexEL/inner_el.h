/// @file inner_el.h

#ifndef INNER_EL_H
#define INNER_EL_H

// #include <Rcpp.h>
#include <RcppEigen.h>
// #include "block_outer.h" // columnwise outer product
#include "adj_G.h" // for support adjustment
#include "log_star.h"

// [[Rcpp::depends(RcppEigen)]]

namespace flexEL {

  using namespace Eigen;

  /// Maximum relative error between two vectors.
  ///
  /// Calculates
  /// ```
  /// max |x1 - x2| / (|x1 + x2| + .1)
  /// ```
  /// The constant in the denominator is for situations in which some of thex's are very small in absolute value.
  ///
  /// @param x1[in] First vector.
  /// @param x2[in] Second vector.
  /// @return The maximum relative error.
  inline double max_rel_err(const Ref<const VectorXd>& x1,
			    const Ref<const VectorXd>&x2) {
    return ((x1 - x2).array().abs() /
    	    ((x1 + x2).array().abs() + 0.1)).maxCoeff();
  }
  
  /// Empirical likelihood for the general moment specification.
  class GenEL {
    
  private:
    
    // dim of G matrix
    int n_obs_; // number of columns
    int n_eqs_; // number of rows    
    // constants for logstar
    double trunc_, aa_, bb_, cc_; 
    // internal variables for lambda_nr
    MatrixXd Gaug_; // augmented G matrix with space for support adjustment 
    VectorXd lambda_new_; // new lambda in Newton-Raphson iterations
    MatrixXd GGt_; //  G times G transpose
    VectorXd Glambda_; // G times lambda
    // ArrayXd Gl11_; // values of omegas before standardization
    VectorXd Q1_; // 1st derivative of the objective function
    MatrixXd Q2_; // 2nd derivative of the objective function
    LLT<MatrixXd> Q2llt_; // robust Cholesky decomposition of Q2
    VectorXd rho_; // placeholder in LambdaNR
    // internal variables for logel
    VectorXd lambda_; // internal lambda for Newton-Raphson iterations
    VectorXd omega_; //  internal weights for empirical distribution
    // MatrixXd GGt_used_;
    // VectorXd rel_err_; //relative error
    
    // internals for Newton-Raphson options
    int max_iter_; // maximum number of iterations
    int nr_iter_; // actual number of iterations
    double rel_tol_; // relative tolerance
    double nr_err_; // actual relative error
    bool supp_adj_; // whether do support adjustment or not
    double supp_a_; // tuning parameter for support currection
    VectorXd lambda0_; // initial lambda in Newton-Raphson iterations

    // for support modification
    int n_obs1_; // n_obs1_ = n_obs_+1: for memory allocation
    int n_obs2_; // n_obs2_ = n_obs_+support: used n_obs_ in calculation
        
    // log_star and its derivatives for the EL dual problem
    double log_star(double x) {return ::flexEL::log_star(x, trunc_, aa_, bb_, cc_);}
    double log_star1(double x) {return ::flexEL::log_star1(x, trunc_, aa_, bb_);}
    double log_star2(double x) {return ::flexEL::log_star2(x, trunc_, aa_);}

    // implementation versions of functions to avoid defining pointer to matrix
    void omega_hat_impl(Ref<VectorXd> omega,
			const Ref<const MatrixXd>& Geff,
			const Ref<const VectorXd>& lambda);
    void lambda_nr_impl(Ref<VectorXd> lambda,
			const Ref<const MatrixXd>& Geff);
      
  public:

    /// Constructor.
    GenEL(int n_obs, int n_eqs);

    /// Set the maximum number of Newton-Raphson iterations.
    void set_max_iter(int max_iter);
    /// Set the relative tolerance for the Newton-Raphson algorithm.
    void set_rel_tol(double rel_tol);
    /// Set the support adjustment flag.
    void set_supp_adj(bool supp_adj, double a);
    void set_supp_adj(bool supp_adj);
    /// Set the initial value for Newton-Raphson algorithm.
    void set_lambda0(const Ref<const VectorXd>& lambda0);
    /// Get the diagnostics for the last Newton-Raphson run.
    void get_diag(int& nr_iter, double& nr_err);

    /// Solve the dual problem via Newton-Raphson algorithm.
    void lambda_nr(Ref<VectorXd> lambda, const Ref<const MatrixXd>& G);
    /// Calculate the profile probability weights.
    void omega_hat(Ref<VectorXd> omega, const Ref<const MatrixXd>& G,
		   const Ref<const VectorXd>& lambda);
    /// Calculate the empirical loglikelihood given the probability weights.
    double logel_omega(const Ref<const VectorXd>& omega);
    /// Calculate the empirical loglikelihood.
    double logel(const Ref<const MatrixXd>& G);
    /// Calculate the gradient of the empirical loglikelihood.
    void logel_grad(Ref<MatrixXd> dldG, const Ref<const MatrixXd>& G,
		    const Ref<const VectorXd>& omega,
		    const Ref<const VectorXd>& lambda);
    double logel_grad(Ref<MatrixXd> dldG, const Ref<const MatrixXd>& G);
  };


  /// @param[in] n_obs    Number of observations.
  /// @param[in] n_eqs    Number of estimating equations.
  inline GenEL::GenEL(int n_obs, int n_eqs) {
    // assign internal values
    n_obs_ = n_obs;
    n_eqs_ = n_eqs;
    // initialization of log_star constants
    trunc_ = 1.0 / n_obs_;
    aa_ = -.5 * n_obs_*n_obs_;
    bb_ = 2.0 * n_obs_;
    cc_ = -1.5 - log(n_obs_);
    n_obs1_ = n_obs_+1; // space for augmented G
    // support adjustment
    supp_adj_ = false;
    n_obs2_ = n_obs_+supp_adj_;
    // memory allocation
    omega_ = VectorXd::Constant(n_obs1_, 1.0/(double)n_obs1_); // Initialize to 1/n_obs_
    Gaug_ = MatrixXd::Zero(n_eqs_,n_obs1_); // augmented G
    // GGt_ = MatrixXd::Zero(n_eqs_,n_obs1_*n_eqs_);
    GGt_= MatrixXd::Zero(n_eqs_, n_eqs_); // columnwise outer product
    lambda_ = VectorXd::Zero(n_eqs_); // Initialize to all 0's
    lambda_new_ = VectorXd::Zero(n_eqs_);
    Q1_ = VectorXd::Zero(n_eqs_);
    Q2_ = MatrixXd::Zero(n_eqs_,n_eqs_);
    Glambda_ = VectorXd::Zero(n_obs1_);
    // Gl11_ = ArrayXd::Zero(n_obs1_);
    rho_ = VectorXd::Zero(n_obs1_);
    Q2llt_.compute(MatrixXd::Identity(n_eqs_,n_eqs_));
    // initialize NR options
    nr_iter_ = 0;
    nr_err_ = 0.0;
    set_max_iter(100);
    set_rel_tol(1e-7);
    set_supp_adj(false);
    set_lambda0(VectorXd::Zero(n_eqs_));
  }

  /// @param[in] max_iter Maximum number of Newton-Raphson iterations.
  inline void GenEL::set_max_iter(int max_iter) {
    max_iter_ = max_iter;
    return;
  }

  /// @param[in] rel_tol Relative tolerance for the Newton-Raphson algorithm.
  inline void GenEL::set_rel_tol(double rel_tol) {
    rel_tol_ = rel_tol;
    return;
  }

  /// Reference: J. Chen, A. M. Variyath, and B. Abraham. Adjusted empirical likelihood and its properties. Journal of Computational and Graphical Statistics, 17(2):426–443, 2008).
  ///
  /// @param[in] supp_adj Whether or not to enable support adjustment.
  /// @param[in] a Support adjustment factor.  Defaults to `max(1.0, log(n_obs)/2)`.
  inline void GenEL::set_supp_adj(bool supp_adj, double a) {
    supp_adj_ = supp_adj;
    supp_a_ = a;
    n_obs2_ = n_obs_+supp_adj_;
    return;
  }
  inline void GenEL::set_supp_adj(bool supp_adj) {
    set_supp_adj(supp_adj, std::max(1.0,0.5*log(n_obs_)));
    return;
  }

  /// @param[in] lambda0 Initialization vector of size `n_eqs`.
  inline void GenEL::set_lambda0(const Ref<const VectorXd>& lambda0) {
    lambda0_ = lambda0;
    return;
  }

  /// @param[out] nr_iter Number of Newton-Raphson iterations.
  /// @param[out] nr_err Maximum relative difference between elements of `lambda` in the last two Newton-Raphson steps.
  inline void GenEL::get_diag(int& nr_iter, double& nr_err) {
    nr_iter = nr_iter_;
    nr_err = nr_err_;
    return;
  }

  inline void GenEL::lambda_nr_impl(Ref<VectorXd> lambda,
				    const Ref<const MatrixXd>& Geff) {
  
    lambda = lambda0_; // set to initial value
    // newton-raphson loop
    int ii, jj;
    for(ii=0; ii<max_iter_; ii++) {
      // Q1 and Q2
      Glambda_.noalias() = lambda.transpose() * Geff;
      Glambda_= 1.0 - Glambda_.array();
      // std::cout << "Glambda = \n" << Glambda_ << std::endl;
      Q2_.fill(0.0);
      for(jj=0; jj<n_obs2_; jj++) {
	rho_(jj) = log_star1(Glambda_(jj));
	GGt_.noalias() = Geff.col(jj) * Geff.col(jj).transpose();
	Q2_ -= log_star2(Glambda_(jj)) * GGt_;
	// Q2_ += log_star2(Glambda_(jj)) * GGt_.block(0,jj*n_eqs_,n_eqs_,n_eqs_);
      }
      // std::cout << "Q2 = \n" << Q2_ << std::endl;
      // update lambda
      Q1_.noalias() = Geff * rho_.head(n_obs2_);
      // std::cout << "Q1 = \n" << Q1_ << std::endl;      
      Q2llt_.compute(Q2_);
      Q2llt_.solveInPlace(Q1_);
      // std::cout << "solve(Q2, Q1) = \n" << Q1_ << std::endl;
      lambda_new_ = lambda - Q1_;
      // std::cout << "lambda_new = \n" << lambda_new_ << std::endl;
      // lambda_new_.noalias() = lambda_old_ - Q2llt_.solve(Q1_);
      nr_err_ = max_rel_err(lambda, lambda_new_);
      lambda = lambda_new_; // complete cycle
      if (nr_err_ < rel_tol_) {
	break;
      }
    }
    nr_iter_ = ii;
    return;
  }


  /// @param[out] lambda Vector of length `n_eqs` containing the candidate solution to the dual optimization problem.
  /// @param[in] G Moment matrix of size `n_eqs x n_obs`.
  ///
  /// @note The internal version of this fnction `lambda_nr_impl` uses the "effective" `G` which has `n_obs + supp_adj` columns.  This is a workaround to assigning a pointer to an Eigen matrix which is specified by a conditional (`supp_adj`).
  inline void GenEL::lambda_nr(Ref<VectorXd> lambda,
			       const Ref<const MatrixXd>& G) {
    if(supp_adj_) {
      Gaug_.leftCols(n_obs_) = G;
      adj_G(Gaug_, supp_a_);
      lambda_nr_impl(lambda, Gaug_);
    } else {
      lambda_nr_impl(lambda, G);
    }
  }

  inline void GenEL::omega_hat_impl(Ref<VectorXd> omega,
				    const Ref<const MatrixXd>& Geff,
				    const Ref<const VectorXd>& lambda) {
    // std::cout << "Geff = \n" << Geff << std::endl;
    Glambda_.head(n_obs2_).noalias() = lambda.transpose() * Geff;
    omega = (1.0-Glambda_.head(n_obs2_).array()).inverse();
    // note: when lambda is the optimal value denominator is n_obs...
    omega.array() /= omega.array().sum();
    return;
  }

  
  /// @param[out] omega Probability vector of length `n_obs + supp_adj`.
  /// @param[in] G Moment matrix of size `n_eqs x n_obs` _or_ of size `n_eqs x (n_obs + 1)`.  When `supp_adj && G.cols() == n_obs`, performs the support adjustment.  Otherwise, assumes that `G` has already been support corrected./// @param[in] lambda Dual problem vector of size `n_eqs`.  
  /// @note The internal version of this fnction `omega_hat_impl` uses the "effective" `G` which has `n_obs + supp_adj` columns.  This is a workaround to assigning a pointer to an Eigen matrix which is specified by a conditional (`supp_adj`).
  inline void GenEL::omega_hat(Ref<VectorXd> omega,
			       const Ref<const MatrixXd>& G,
			       const Ref<const VectorXd>& lambda) {
    if(supp_adj_ && (G.cols() == n_obs_)) {
      Gaug_.leftCols(n_obs_) = G;
      adj_G(Gaug_, supp_a_);
      omega_hat_impl(omega, Gaug_, lambda);
    } else {
      omega_hat_impl(omega, G, lambda);
    }
    // Ref<MatrixXd> G_ = (supp_adj_ && (G.cols() == n_obs_)) ? Gaug_ : G;
    // Glambda_.head(n_obs2_).noalias() = lambda.transpose() * G_;
    // omega = (1.0-Glambda_.head(n_obs2_).array()).inverse();
    // // note: when lambda is the optimal value denominator is n_obs...
    // omega.array() /= omega.array().sum();
    return;
  }

  /// @param[in] omega Probability vector of length `n_obs + supp_adj`.
  /// @return Value of empirical loglikelihood, which is `sum(log(omega))`.
  inline double GenEL::logel_omega(const Ref<const VectorXd>& omega) {
    return omega.array().log().sum();
  }

  /// @param[in] G Moment matrix of size `n_eqs x n_obs`.
  /// @return Value of empirical loglikelihood, which is `sum(log(omega))`.
  inline double GenEL::logel(const Ref<const MatrixXd>& G) {
    lambda_nr(lambda_, G); // newton-raphson step
    if(supp_adj_) {
      omega_hat(omega_, Gaug_, lambda_);
    } else {
      omega_hat(omega_.head(n_obs2_), G, lambda_);
    }
    return logel_omega(omega_.head(n_obs2_));
  }

  /// @param[out] dldG Matrix of size `n_eqs x n_obs` containing the gradient of logel with respect to G.
  /// @param[in] G Moment matrix of size `n_eqs x n_obs`.
  /// @param[out] omega Probability vector of length `n_obs + supp_adj`.
  /// @param[in] lambda Dual problem vector of size `n_eqs`.
  inline void GenEL::logel_grad(Ref<MatrixXd> dldG,
				const Ref<const MatrixXd>& G,
				const Ref<const VectorXd>& omega,
				const Ref<const VectorXd>& lambda) {
    dldG.noalias() = lambda * omega.head(n_obs_).transpose();
    if(supp_adj_) {
      dldG.colwise() -= omega(n_obs1_-1) * supp_a_/n_obs_ * lambda;
    }
    dldG *= n_obs2_;
    return;
  }

  /// @param[out] dldG Matrix of size `n_eqs x n_obs` containing the gradient of logel with respect to G.
  /// @param[in] G Moment matrix of size `n_eqs x n_obs`.
  /// @return For convenience, returns the value of the empirical loglikelihood.
  inline double GenEL::logel_grad(Ref<MatrixXd> dldG,
				  const Ref<const MatrixXd>& G) {
    double ll = logel(G);
    // use internal values of omega_ and lambda_
    logel_grad(dldG, G, omega_, lambda_);
    return ll;
  }

} // end namespace flexEL

#endif

// --- scratch -----------------------------------------------------------------

// // set options functions

// /**
//  * @brief Set tolerance values for NR, support adjustment option and tuning parameter, and initial value of lambda.
//  * 
//  * @param[in] max_iter    Maximum number of iterations.
//  * @param[in] rel_tol     Relative tolerance to control the convergence of lambda.
//  * @param[in] supp_adj        Whether to conduct support adjustment.
//  * @param[in] supp_a      Tuning parameter for support adjustment (see J. Chen, A. M. Variyath, and B. Abraham. Adjusted empirical likelihood and its properties. Journal of Computational and Graphical Statistics, 17(2):426–443, 2008).
//  * @param[in] lambda0     Initial value for lambda.
//  */
// inline void flexEL::InnerEL::set_opts(const int& max_iter, const double& rel_tol, 
//                                       const bool& supp_adj, const double& supp_a, 
//                                       const Ref<const VectorXd>& lambda0) {
//   max_iter_ = max_iter;
//   rel_tol_ = rel_tol;
//   supp_adj_ = supp_adj;
//   supp_a_ = supp_a;
//   n_obs2_ = n_obs_+supp_adj_;
//   lambda0_ = lambda0;
// }

// /**
//  * @brief Set tolerance values for NR, support adjustment option, and initial value of lambda.
//  * 
//  * @param[in] max_iter    Maximum number of iterations.
//  * @param[in] rel_tol     Relative tolerance to control the convergence of lambda.
//  * @param[in] supp_adj        Whether to conduct support adjustment.
//  * @param[in] lambda0     Initial value for lambda.
//  */
// inline void flexEL::InnerEL::set_opts(const int& max_iter, const double& rel_tol, 
//                                       const bool& supp_adj, 
//                                       const Ref<const VectorXd>& lambda0) {
//   max_iter_ = max_iter;
//   rel_tol_ = rel_tol;
//   supp_adj_ = supp_adj;
//   supp_a_ = std::max(1.0,0.5*log(n_obs_));
//   n_obs2_ = n_obs_+supp_adj_;
//   lambda0_ = lambda0;
// }

// /**
//  * @brief Set tolerance values for NR, support adjustment option. Initial value of lambda is omitted and is default to be a vector of zeros.
//  * 
//  * @param[in] max_iter    Maximum number of iterations.
//  * @param[in] rel_tol     Relative tolerance to control the convergence of lambda.
//  * @param[in] supp_adj        Whether to conduct support adjustment.
//  * @param[in] supp_a      Tuning parameter for support adjustment (see J. Chen, A. M. Variyath, and B. Abraham. Adjusted empirical likelihood and its properties. Journal of Computational and Graphical Statistics, 17(2):426–443, 2008).
//  */
// inline void flexEL::InnerEL::set_opts(const int& max_iter, const double& rel_tol, 
//                                       const bool& supp_adj, const double& supp_a) {
//   max_iter_ = max_iter;
//   rel_tol_ = rel_tol;
//   supp_adj_ = supp_adj;
//   supp_a_ = supp_a;
//   n_obs2_ = n_obs_+supp_adj_;
// }

// /**
//  * @brief Set tolerance values for NR, support adjustment option. Initial value of lambda is omitted and is default to be a vector of zeros.
//  * 
//  * @param[in] max_iter    Maximum number of iterations.
//  * @param[in] rel_tol     Relative tolerance to control the convergence of lambda.
//  * @param[in] supp_adj        Whether to conduct support adjustment.
//  */
// inline void flexEL::InnerEL::set_opts(const int& max_iter, const double& rel_tol, 
//                                       const bool& supp_adj) {
//   max_iter_ = max_iter;
//   rel_tol_ = rel_tol;
//   supp_adj_ = supp_adj;
//   supp_a_ = std::max(1.0,0.5*log(n_obs_));
//   n_obs2_ = n_obs_+supp_adj_;
// }

// /**
//  * @brief Set support adjustment option and tuning parameter only. 
//  * 
//  * @param[in] supp_adj      Whether to conduct support adjustment.
//  * @param[in] supp_a    Tuning parameter for support adjustment (see J. Chen, A. M. Variyath, and B. Abraham. Adjusted empirical likelihood and its properties. Journal of Computational and Graphical Statistics, 17(2):426–443, 2008).
//  */
// inline void flexEL::InnerEL::set_opts(const bool& supp_adj, const double& supp_a) {
//   supp_adj_ = supp_adj;
//   supp_a_ = supp_a;
//   n_obs2_ = n_obs_+supp_adj_;
// }

// /**
//  * @brief Set support adjustment option only. 
//  * 
//  * @param[in] supp_adj   Whether to conduct support adjustment.
//  */
// inline void flexEL::InnerEL::set_opts(const bool& supp_adj) {
//   supp_adj_ = supp_adj;
//   supp_a_ = std::max(1.0,0.5*log(n_obs_));
//   n_obs2_ = n_obs_+supp_adj_;
// }

// // maximum relative error in lambda
// // note that a small constant is added to the denominator for situations in which some of the lambda's are converging to a very small number.
// inline double flexEL::GenEL::max_rel_err() {
//   // // TODO: added for numerical stability, what is a good tolerance to use ?
//   // if ((lambda_new_ - lambda_old_).array().abs().maxCoeff() < 1e-10) return(0);
//   // rel_err_ = ((lambda_new_ - lambda_old_).array() / (lambda_new_ + lambda_old_).array()).abs();
//   // return(rel_err_.maxCoeff());
//   return ((lambda_new_ - lambda_old_).array().abs() /
// 	  ((lambda_new_ + lambda_old_).array().abs() + 0.1)).maxCoeff();
// }

// /**
//  * @brief Evaluate omegas based on current G and lambda.
//  */
// inline void flexEL::InnerEL::EvalOmegas() {
//   // G and lambdaNew must have been assigned
//   if (lambda_new_ != lambda_new_) { // if lambdaNew is NaN 
//     for (int ii=0; ii<n_obs2_; ii++) {
//       omegas_(ii) = std::numeric_limits<double>::quiet_NaN();
//     }
//   }
//   else {
//     // if (supp_adj_) adj_G(G_,supp_a_); // calculate based on adjusted G
//     MatrixXd G_used_ = G_.block(0,0,n_eqs_,n_obs2_); // TODO: new allocation right now..
//     Glambda_.head(n_obs2_).noalias() = lambda_new_.transpose() * G_used_;
//     Gl11_.head(n_obs2_) = 1.0/(1.0-Glambda_.head(n_obs2_).array());
//     omegas_.head(n_obs2_).array() = Gl11_.head(n_obs2_).array() / Gl11_.head(n_obs2_).sum(); // in fact, Gl11.sum() should be equal to n_obs
//   }
// }

// /**
//  * @brief Calculate LogEL using omegas.
//  * 
//  * @return log empirical likelihood.
//  */
// inline double flexEL::InnerEL::LogEL() {
//   // if omegas are NaN, return -Inf
//   if (omegas_.head(n_obs2_) != omegas_.head(n_obs2_)) return -INFINITY;
//   else return(omegas_.head(n_obs2_).array().log().sum());
// }

// /**
//  * @brief Calculate LogEL and gradient matrix dldG
//  * 
//  * @param[out] logel   The value of log EL evaluated at the current G.
//  * @param[out] dldG    The first derivative of log EL w.r.t. G evaluated at the current G.
//  */
// inline void flexEL::InnerEL::LogELGrad(double& logel, MatrixXd& dldG) {
//   if (omegas_.head(n_obs2_) != omegas_.head(n_obs2_)) {
//     logel = -INFINITY;
//     dldG = MatrixXd::Zero(n_eqs_,n_obs2_);
//   }
//   logel = omegas_.head(n_obs2_).array().log().sum();
//   if (supp_adj_ == false) {
//     dldG = n_obs_ * lambda_new_ * omegas_.head(n_obs_).transpose();
//   }
//   else {
//     dldG = n_obs2_ * lambda_new_ * omegas_.head(n_obs_).transpose() - 
//       n_obs2_ * omegas_(n_obs2_-1) * supp_a_/n_obs_ * lambda_new_ * VectorXd::Ones(n_obs_).transpose();
//   }
// }

// // setters

// /**
//  * @brief Set the value of lambda (e.g. to be used directly to calculate omegas).
//  * 
//  * @param[in] lambda   A numeric vector.
//  */
// inline void flexEL::InnerEL::set_lambda(const Ref<const VectorXd>& lambda) {
//   // lambda_old_ = lambda;
//   lambda_new_ = lambda;
// }

// /**
//  * @brief Set the value of omegas (e.g. to be used directly to calculate log EL).
//  * 
//  * @param[in] omegas   A numeric probability vector (which sums to 1).
//  */
// inline void flexEL::InnerEL::set_omegas(const Ref<const VectorXd>& omegas) {
//   omegas_.head(n_obs2_) = omegas; 
// }

// /**
//  * @brief Set the value of G (e.g. to be used directly to calculate lambda or log EL).
//  * 
//  * @param[in] G   A numeric matrix of dimension n_eqs x n_obs.
//  */
// inline void flexEL::InnerEL::set_G(const Ref<const MatrixXd>& G) {
//   G_.block(0,0,n_eqs_,n_obs_) = G;
//   if (supp_adj_) adj_G(G_,supp_a_);
// }

// // getters

// /**
//  * @brief Get the value of lambda.
//  */
// inline VectorXd flexEL::InnerEL::get_lambda() {
//   return(lambda_new_);
// }

// /**
//  * @brief Get the value of omegas.
//  */
// inline VectorXd flexEL::InnerEL::get_omegas() {
//   return(omegas_.head(n_obs2_));
// }

// /**
//  * @brief Get the value of G.
//  */
// inline MatrixXd flexEL::InnerEL::get_G() {
//   return(G_.block(0,0,n_eqs_,n_obs2_));
// }

// /**
//  * @brief Get the reference of original dimension G (not including the adjusted row if exists).
//  */
// inline Ref<MatrixXd> flexEL::InnerEL::get_ref_G() {
//   return Ref<MatrixXd>(G_.block(0,0,n_eqs_,n_obs2_-supp_adj_));
// }
