/// @file gen_el.h

#ifndef FLEXEL_GEN_EL_H
#define FLEXEL_GEN_EL_H

#include "utils.h"
#include <Eigen/Dense>

namespace flexEL {

  using namespace Eigen;
  
  /// Empirical likelihood for the general moment specification.
  template <class Type>
  class GenEL {
    
  private:
    
    // dim of G matrix
    int n_obs_; // number of columns
    int n_eqs_; // number of rows    
    // for support modification
    int n_obs1_; // n_obs1_ = n_obs_+1: for memory allocation
    int n_obs2_; // n_obs2_ = n_obs_+support: used n_obs_ in calculation
    // logstar parameters
    Vector_t<Type> tstar_;
    Vector_t<Type> astar_;
    Vector_t<Type> bstar_;
    // Vector_t<Type> cstar_;
    Vector_t<Type> norm_weights_;
    // Type trunc_, aa_, bb_, cc_; 
    // internal variables for Newton-Raphson 
    Matrix_t<Type> G_aug_; // augmented G matrix with space for support adjustment 
    Vector_t<Type> lambda_new_; // new lambda in Newton-Raphson iterations
    Matrix_t<Type> GGt_; //  G times G transpose
    Vector_t<Type> Glambda_; // G times lambda
    // ArrayXd Gl11_; // values of omegas before standardization
    Vector_t<Type> Q1_; // 1st derivative of the objective function
    Matrix_t<Type> Q2_; // 2nd derivative of the objective function
    LLT<Matrix_t<Type> > Q2llt_; // robust Cholesky decomposition of Q2
    Vector_t<Type> rho_; // placeholder in LambdaNR
    // internal variables for logel
    Vector_t<Type> lambda_; // internal lambda for Newton-Raphson iterations
    Vector_t<Type> omega_; //  internal weights for empirical distribution
    // Matrix_t<Type> GGt_used_;
    // Vector_t<Type> rel_err_; //relative error
    
    // internals for Newton-Raphson options
    int max_iter_; // maximum number of iterations
    int nr_iter_; // actual number of iterations
    Type rel_tol_; // relative tolerance
    Type nr_err_; // actual relative error
    bool supp_adj_; // whether do support correction or not
    Type supp_a_; // tuning parameter for support correction
    Type weight_adj_; // weight of the additional observation under support correction
    Vector_t<Type> lambda0_; // initial lambda in Newton-Raphson iterations

    // implementation versions of functions to avoid defining pointer to matrix
    void omega_hat_impl(RefVector_t<Type> omega,
                        cRefVector_t<Type>& lambda,
                        cRefMatrix_t<Type>& G,
                        cRefVector_t<Type>& norm_weights);
    void lambda_nr_impl(RefVector_t<Type> lambda,
                        cRefMatrix_t<Type>& G,
                        cRefVector_t<Type>& norm_weights);
    Type logel_omega_impl(cRefVector_t<Type>& omega,
			  cRefVector_t<Type>& norm_weights,
			  Type sum_weights);
    Type logel_full_impl(RefVector_t<Type> omega,
			 RefVector_t<Type> lambda,
			 cRefMatrix_t<Type>& G,
			 cRefVector_t<Type>& norm_weights,
			 Type sum_weights);
    Type logel_impl(cRefMatrix_t<Type>& G,
		    cRefVector_t<Type>& norm_weights,
		    Type sum_weights);
      
  public:

    /// Constructor.
    GenEL(int n_obs, int n_eqs);

    /// Set the maximum number of Newton-Raphson iterations.
    void set_max_iter(int max_iter);
    /// Set the relative tolerance for the Newton-Raphson algorithm.
    void set_rel_tol(Type rel_tol);
    /// Set the support adjustment flag and support adjustment parameters.
    void set_supp_adj(bool supp_adj, Type a, Type weight_adj);
    /// Set the support adjustment flag.
    void set_supp_adj(bool supp_adj);
    /// Set the support adjustment tuning parameter.
    void set_supp_adj_a(Type a);
    /// Set the weight for the additional fake observation under support adjustment.
    void set_weight_adj(Type weight_adj);
    /// Set the initial value for Newton-Raphson algorithm.
    void set_lambda0(cRefVector_t<Type>& lambda0);
    /// Get the diagnostics for the last Newton-Raphson run.
    void get_diag(int& nr_iter, Type& nr_err);
    /// Get support adjustment flag.
    bool get_supp_adj();
    /// Get number of observations.
    int get_n_obs();
    /// Get number of estimating equations.
    int get_n_eqs();
    /// Get the maximum number of iterations.
    int get_max_iter();
    /// Get the relative tolerance.
    Type get_rel_tol();

    /// Moment matrix with support adjustment.
    Ref<const Matrix_t<Type> > supp_G(cRefMatrix_t<Type>& G);
    /// Normalized weights with support adjustment.
    Ref<const Vector_t<Type> > supp_norm_weights(cRefVector_t<Type>& weights);
    /// Normalized weights and total weight with support adjustment.
    Ref<const Vector_t<Type> > supp_norm_weights(cRefVector_t<Type>& weights, Type& sum_weights);

    
    /// Solve the dual problem via Newton-Raphson algorithm.
    void lambda_nr(RefVector_t<Type> lambda, 
                   cRefMatrix_t<Type>& G,
                   cRefVector_t<Type>& weights);
    // cRefVector_t<Type>& norm_weights);
    /// Check whether Newton-Raphson algorithm has converged.
    bool has_converged_nr();
    /// Calculate the profile probability weights.
    void omega_hat(RefVector_t<Type> omega,
                   cRefVector_t<Type>& lambda,
                   cRefMatrix_t<Type>& G,
                   cRefVector_t<Type>& weights);
    // cRefVector_t<Type>& norm_weights);
    /// Calculate the empirical loglikelihood given the probability weights.
    Type logel_omega(cRefVector_t<Type>& omega,
		     cRefVector_t<Type>& weights);
    /// Calculate the empirical loglikelihood returning intermediate computations.
    Type logel_full(RefVector_t<Type> omega,
		    RefVector_t<Type> lambda,
		    cRefMatrix_t<Type>& G,
		    cRefVector_t<Type>& weights);
    /// Calculate the unweighted empirical loglikelihood.
    Type logel(cRefMatrix_t<Type>& G);
    /// Calculate the weighted empirical loglikelihood.
    Type logel(cRefMatrix_t<Type>& G,
	       cRefVector_t<Type>& weights);
    /// Calculate the gradient of the unweighted empirical loglikelihood.
    Type logel_grad(RefMatrix_t<Type> dldG, 
		    cRefMatrix_t<Type>& G);
    /// Calculate the gradient of the weighted empirical loglikelihood.
    Type logel_grad(RefMatrix_t<Type> dldG, 
		    cRefMatrix_t<Type>& G,
		    cRefVector_t<Type>& weights);
    /// Calculate the gradient of the empirical loglikelihood.
    void logel_grad(RefMatrix_t<Type> dldG,
                    cRefVector_t<Type>& omega,
                    cRefVector_t<Type>& lambda,
                    Type sum_weights);
  };

  /// @param[in] n_obs    Number of observations.
  /// @param[in] n_eqs    Number of estimating equations.
  template <class Type>
  inline GenEL<Type>::GenEL(int n_obs, int n_eqs) {
    // assign internal values
    n_obs_ = n_obs;
    n_eqs_ = n_eqs;
    n_obs1_ = n_obs_+1; // space for augmented G
    tstar_ = Vector_t<Type>::Zero(n_obs1_);
    astar_ = Vector_t<Type>::Zero(n_obs1_);
    bstar_ = Vector_t<Type>::Zero(n_obs1_);
    norm_weights_ = Vector_t<Type>::Zero(n_obs1_);
    // memory allocation
    omega_ = Vector_t<Type>::Constant(n_obs1_, 1.0/Type(n_obs1_)); // Initialize to 1/n_obs_
    G_aug_ = Matrix_t<Type>::Zero(n_eqs_,n_obs1_); // augmented G
    GGt_= Matrix_t<Type>::Zero(n_eqs_, n_eqs_); // columnwise outer product
    lambda_ = Vector_t<Type>::Zero(n_eqs_); // Initialize to all 0's
    lambda_new_ = Vector_t<Type>::Zero(n_eqs_);
    Q1_ = Vector_t<Type>::Zero(n_eqs_);
    Q2_ = Matrix_t<Type>::Zero(n_eqs_,n_eqs_);
    Glambda_ = Vector_t<Type>::Zero(n_obs1_);
    rho_ = Vector_t<Type>::Zero(n_obs1_);
    Q2llt_.compute(Matrix_t<Type>::Identity(n_eqs_,n_eqs_));
    // initialize NR options
    nr_iter_ = 0;
    nr_err_ = 0.0;
    set_max_iter(100);
    set_rel_tol(1e-7);
    set_supp_adj(false);
    set_weight_adj(1.0);
    set_lambda0(Vector_t<Type>::Zero(n_eqs_));
  }

  /// @param[in] max_iter Maximum number of Newton-Raphson iterations.
  template <class Type>
  inline void GenEL<Type>::set_max_iter(int max_iter) {
    max_iter_ = max_iter;
    return;
  }

  /// @param[in] rel_tol Relative tolerance for the Newton-Raphson algorithm.
  template <class Type>
  inline void GenEL<Type>::set_rel_tol(Type rel_tol) {
    rel_tol_ = rel_tol;
    return;
  }

  /// Reference: J. Chen, A. M. Variyath, and B. Abraham. Adjusted empirical likelihood and its properties. Journal of Computational and Graphical Statistics, 17(2):426–443, 2008).
  ///
  /// @param[in] supp_adj   Whether or not to enable support adjustment.
  /// @param[in] a          Support adjustment factor.  Defaults to `max(1.0, log(n_obs)/2)`.
  /// @param[in] weight_adj Weight for the additional observation when calculating weighted log EL.
  template <class Type>
  inline void GenEL<Type>::set_supp_adj(bool supp_adj, Type a, Type weight_adj) {
    supp_adj_ = supp_adj;
    supp_a_ = a;
    weight_adj_ = weight_adj;
    n_obs2_ = n_obs_+supp_adj_;
    return;
  }
  
  /// Reference: J. Chen, A. M. Variyath, and B. Abraham. Adjusted empirical likelihood and its properties. Journal of Computational and Graphical Statistics, 17(2):426–443, 2008).
  ///
  /// @param[in] supp_adj   Whether or not to enable support adjustment.
  template <class Type>
  inline void GenEL<Type>::set_supp_adj(bool supp_adj) {
    set_supp_adj(supp_adj, std::max(1.0,0.5*log(n_obs_)), 1.0);
    return;
  }
  
  /// Reference: J. Chen, A. M. Variyath, and B. Abraham. Adjusted empirical likelihood and its properties. Journal of Computational and Graphical Statistics, 17(2):426–443, 2008).
  ///
  /// @param[in] a          Support adjustment factor.  Defaults to `max(1.0, log(n_obs)/2)`.
  template <class Type>
  inline void GenEL<Type>::set_supp_adj_a(Type a) {
    supp_a_ = a;
    return;
  }
  
  /// @param[in] weight_adj Weight for the additional observation when calculating weighted log EL.
  template <class Type>
  inline void GenEL<Type>::set_weight_adj(Type weight_adj) {
    weight_adj_ = weight_adj;
    return;
  }

  /// @param[in] lambda0 Initial value of lambda of length `n_eqs`.
  template <class Type>
  inline void GenEL<Type>::set_lambda0(cRefVector_t<Type>& lambda0) {
    lambda0_ = lambda0;
    return;
  }

  template <class Type>
  inline bool GenEL<Type>::get_supp_adj() {
    return supp_adj_;
  }

  template <class Type>
  inline int GenEL<Type>::get_n_obs() {
    return n_obs_;
  }
  
  template <class Type>
  inline int GenEL<Type>::get_n_eqs() {
    return n_eqs_;
  }
  
  template <class Type>
  inline int GenEL<Type>::get_max_iter() {
    return max_iter_;
  }
  
  template <class Type>
  inline Type GenEL<Type>::get_rel_tol() {
    return rel_tol_;
  }

  /// @param[out] nr_iter Number of Newton-Raphson iterations.
  /// @param[out] nr_err Maximum relative difference between elements of `lambda` in the last two Newton-Raphson steps.
  template <class Type>
  inline void GenEL<Type>::get_diag(int& nr_iter, Type& nr_err) {
    nr_iter = nr_iter_;
    nr_err = nr_err_;
    return;
  }

  /// @return Whether the desired error tolerance `rel_tol` has been reached in less than `max_iter` Newton-Raphson steps.
  template <class Type>
  inline bool GenEL<Type>::has_converged_nr() {
    return !((nr_err_ > rel_tol_) && (nr_iter_ == max_iter_));
  }

  /// @param[in] G Moment matrix of size `n_eqs x n_obs` or `n_eqs x (n_obs + 1)`. If `supp_adj = false`, the former is required.  If `supp_adj = true` and the former is provided, support adjustment is performed.  If `supp_adj = true` and `G.cols() == n_obs + 1`, assumes that support has already been corrected.
  /// return A reference to `G` if no support adjustment is necessary (not needed or already performed); otherwise a reference to `G_aug_`. 
  template <class Type>
  inline Ref<const Matrix_t<Type> > GenEL<Type>::supp_G(cRefMatrix_t<Type>& G) {
    if(supp_adj_ && G.cols() == n_obs_) {
      G_aug_.leftCols(n_obs_) = G;
      adj_G<Type>(G_aug_, supp_a_);
      return Ref<const Matrix_t<Type> >(G_aug_);
    } else {
      return Ref<const Matrix_t<Type> >(G);
    }
  }

  /// @param[in] weights Vector of weights of length `n_obs` or `n_obs+1`.  If `supp_adj = false`, the former is required.  If `supp_adj = true` and the former is provided, the internal value of `weight_adj` is appended.  If `supp_adj = false` and `weights.size() == n_obs+1`, assumes that the last entry of `weight` already contains a weight adjustment.
  /// @param[out] sum_weights Sum of `weights` after potential support adjustment.
  /// @return A reference to a vector of length `n_obs+supp_adj` consisting of the normalized version of the possibly modified `weights`, i.e., entries sum to one.
  template <class Type>
  inline Ref<const Vector_t<Type> > GenEL<Type>::supp_norm_weights(cRefVector_t<Type>& weights, Type& sum_weights) {
    if(supp_adj_ && weights.size() == n_obs_) {
      norm_weights_.head(n_obs_) = weights;
      norm_weights_(n_obs_) = weight_adj_;
      // std::cout << "norm_weights_ before normalization = " << norm_weights_.transpose() << std::endl;
      sum_weights = norm_weights_.sum();
      norm_weights_ /= sum_weights;
      return Ref<const Vector_t<Type> >(norm_weights_);
    } else {
      // std::cout << "supp_norm_weights: only normalize." << std::endl;
      // std::cout << "weights = " << weights.transpose() << std::endl;
      int n_obs2 = weights.size();
      sum_weights = weights.sum();
      norm_weights_.head(n_obs2) = weights;
      norm_weights_.head(n_obs2) /= sum_weights;
      return Ref<const Vector_t<Type> >(norm_weights_.head(n_obs2));
    }
  }

  template <class Type>
  inline Ref<const Vector_t<Type> > GenEL<Type>::supp_norm_weights(cRefVector_t<Type>& weights) {
    Type sum_weights;
    return supp_norm_weights(weights, sum_weights);
  }
  

  /// @param[in/out] lambda Vector of length `n_eqs` containing the initial and final solution to the dual optimization problem.
  /// @param[in] G Moment matrix of size `n_eqs x n_obs` or `n_eqs x (n_obs+1)`.
  /// @param[in] norm_weights Vector of size `G.cols()` of normalized weights (summing to one).
  template <class Type>
  inline void GenEL<Type>::lambda_nr_impl(RefVector_t<Type> lambda,
					  cRefMatrix_t<Type>& G,
					  cRefVector_t<Type>& norm_weights) {
    // printf("Entering lambda_nr_impl()\n");
    // std::cout << "lambda:" << std::endl << lambda.transpose() << std::endl;
    // std::cout << "G:" << std::endl << G << std::endl;
    // std::cout << "norm_weights:" << std::endl << norm_weights.transpose() << std::endl;
    int n_obs2 = G.cols();
    // logstar constants
    // assume the weights are normalized to sum to one
    bstar_.head(n_obs2) = norm_weights.array().inverse();
    astar_.head(n_obs2) = -.5 * bstar_.head(n_obs2).array().square();
    bstar_.head(n_obs2) *= 2.0;
    // newton-raphson loop
    int ii, jj;
    for(ii=0; ii<max_iter_; ii++) {
      // Q1 and Q2
      Glambda_.head(n_obs2).noalias() = lambda.transpose() * G;
      Glambda_.head(n_obs2) = 1.0 - Glambda_.head(n_obs2).array();
      Q2_.fill(0.0);
      for(jj=0; jj<n_obs2; jj++) {
        rho_(jj) = log_star1<Type>(Glambda_(jj),
				   norm_weights(jj), astar_(jj), bstar_(jj));
        GGt_.noalias() = G.col(jj) * G.col(jj).transpose();
        Q2_ -= norm_weights(jj) *
          log_star2<Type>(Glambda_(jj), norm_weights(jj), astar_(jj)) * GGt_;
      }
      // update lambda
      rho_.head(n_obs2).array() *= norm_weights.array();
      Q1_.noalias() = G * rho_.head(n_obs2);
      Q2llt_.compute(Q2_);
      Q2llt_.solveInPlace(Q1_);
      lambda_new_ = lambda - Q1_;
      nr_err_ = max_rel_err<Type>(lambda, lambda_new_);
      lambda = lambda_new_; // complete cycle
      // printf("lambda[%i]:\n", ii);
      // std::cout << lambda.transpose() << std::endl;
      if (nr_err_ < rel_tol_) {
        break;
      }
    }
    nr_iter_ = ii;
    return;
  }

  /// @param[out] omega Probability vector of size `n_obs` or `n_obs+1`.
  /// @param[in] lambda Dual problem vector of size `n_eqs`.
  /// @param[in] G Moment matrix of size `n_eqs x omega.size()`.
  /// @param[in] norm_weights Normalized weight vector of size `omega.size()`.
  template <class Type>
  inline void GenEL<Type>::omega_hat_impl(RefVector_t<Type> omega,
					  cRefVector_t<Type>& lambda,
					  cRefMatrix_t<Type>& G,
					  cRefVector_t<Type>& norm_weights) {
    int n_obs2 = omega.size();
    Glambda_.head(n_obs2).noalias() = lambda.transpose() * G;
    omega = (1.0-Glambda_.head(n_obs2).array()).inverse() * norm_weights.array();
    // note: when lambda is the optimal value denominator is n_obs.
    // however things get more mixed up with support correction,
    // so explicitly normalize for safety.
    omega.array() /= omega.array().sum();
    return;
  }

  /// @param[in] omega Probability vector of length `n_obs + supp_adj`.
  /// @param[in] norm_weights Normalized weight vector of size `omega.size()`.
  /// @return Value of empirical loglikelihood, which is `sum(log(omega))`.
  template <class Type>
  inline Type GenEL<Type>::logel_omega_impl(cRefVector_t<Type>& omega,
					    cRefVector_t<Type>& norm_weights,
					    Type sum_weights) {
    // Type sum_weights;
    // Ref<const Vector_t<Type> > norm_weights_eff = supp_norm_weights(weights, sum_weights);
    // std::cout << "sum_weights = " << sum_weights << std::endl;
    // std::cout << "omega = " << omega.transpose() << std::endl;
    // std::cout << "norm_weights_eff = " << norm_weights_eff.transpose() << std::endl;
    return sum_weights * (omega.array().log() * norm_weights.array()).sum();
  }


  /// @param[out] omega Probability vector of size `n_obs` or `n_obs+1`.
  /// @param[out] lambda Dual problem vector of size `n_eqs`.
  /// @param[in] G Moment matrix of size `n_eqs x omega.size()`.
  /// @param[in] norm_weights Normalized weight vector of size `omega.size()`.
  /// @param[in] sum_weights Sum of weights prior to normalization.
  ///
  /// @return Value of empirical loglikelihood.
  template <class Type>
  inline Type GenEL<Type>::logel_full_impl(RefVector_t<Type> omega,
					   RefVector_t<Type> lambda,
					   cRefMatrix_t<Type>& G,
					   cRefVector_t<Type>& norm_weights,
					   Type sum_weights) {
    lambda = lambda0_; // initial value
    lambda_nr_impl(lambda, G, norm_weights);
    omega_hat_impl(omega, lambda, G, norm_weights);
    return logel_omega_impl(omega, norm_weights, sum_weights);
  }

  /// @param[in] G Moment matrix of size `n_eqs x n_obs` or `n_eqs x (n_obs+1)`.
  /// @param[in] norm_weights Normalized weight vector of size `G.cols()`.
  /// @param[in] sum_weights Sum of weights prior to normalization.
  ///
  /// @return Value of empirical loglikelihood.
  template <class Type>
  inline Type GenEL<Type>::logel_impl(cRefMatrix_t<Type>& G,
				      cRefVector_t<Type>& norm_weights,
				      Type sum_weights) {
    return logel_full_impl(omega_.head(n_obs2_), lambda_, G, norm_weights, sum_weights);
    // lambda_ = lambda0_; // initial value
    // lambda_nr_impl(lambda_, G, norm_weights);
    // omega_hat_impl(omega_.head(n_obs2_), lambda_, G, norm_weights);
    // return logel_omega_impl(omega_.head(n_obs2_), norm_weights, sum_weights);
  }


  /// @param[out] lambda Vector of length `n_eqs` containing the initial and final solution to the dual optimization problem.
  /// @param[in] G Moment matrix of size `n_eqs x n_obs` or `n_eqs x (n_obs + supp_adj)`.  If `supp_adj = false`, the former is required.  If `supp_adj = true` and the former is provided, support adjustment is performed.  If `supp_adj = true` and `G.cols() == n_obs + 1`, assumes that support has already been corrected. 
  /// @param[in] weights Vector of weights of length `n_obs` or `n_obs + 1` following same logic as above.
  ///
  /// @note Due to the difficulty of creating Eigen "pointers", current hack is to explicitly tell lambda_nr_impl whether to use external or internal `G` and `norm_weights`.
  template <class Type>
  inline void GenEL<Type>::lambda_nr(RefVector_t<Type> lambda,
				     cRefMatrix_t<Type>& G,
				     cRefVector_t<Type>& weights) {
    lambda = lambda0_; // initial value
    Ref<const Matrix_t<Type> > G_eff = supp_G(G);
    Ref<const Vector_t<Type> > norm_weights_eff = supp_norm_weights(weights);
    // std::cout << "norm_weights_eff = " << norm_weights_eff.transpose() << std::endl;
    // std::cout << "lambda_in = " << lambda.transpose() << std::endl;
    // std::cout << "G_eff = " << G_eff << std::endl;
    lambda_nr_impl(lambda, G_eff, norm_weights_eff);
    // std::cout << "lambda_out = " << lambda.transpose() << std::endl;    
    return;
  }
  
  /// @param[out] omega Probability vector of length `n_obs + supp_adj`.
  /// @param[in] lambda Dual problem vector of size `n_eqs`.  
  /// @param[in] G Moment matrix of size `n_eqs x n_obs` or `n_eqs x (n_obs + supp_adj)`.  If `supp_adj = false`, the former is required.  If `supp_adj = true` and the former is provided, support adjustment is performed.  If `supp_adj = true` and `G.cols() == n_obs + 1`, assumes that support has already been corrected. 
  /// @param[in] weights Vector of weights of length `n_obs` or `n_obs + supp_adj` following same logic as above.
  template <class Type>
  inline void GenEL<Type>::omega_hat(RefVector_t<Type> omega,
				     cRefVector_t<Type>& lambda,
				     cRefMatrix_t<Type>& G,
				     cRefVector_t<Type>& weights) {
    // cRefVector_t<Type>& norm_weights) {
    Ref<const Matrix_t<Type> > G_eff = supp_G(G);
    Ref<const Vector_t<Type> > norm_weights_eff = supp_norm_weights(weights);
    omega_hat_impl(omega, lambda, G_eff, norm_weights_eff); 
    return;
  }

  /// @param[in] omega Probability vector of length `n_obs + supp_adj`.
  /// @param[in] weights Vector of weights of length `n_obs` or `n_obs + supp_adj`.
  /// @return Value of empirical loglikelihood, which is `sum(log(omega))`.
  template <class Type>
  inline Type GenEL<Type>::logel_omega(cRefVector_t<Type>& omega,
				       cRefVector_t<Type>& weights) {
    Type sum_weights;
    Ref<const Vector_t<Type> > norm_weights_eff = supp_norm_weights(weights, sum_weights);
    return logel_omega_impl(omega, norm_weights_eff, sum_weights);
    // std::cout << "sum_weights = " << sum_weights << std::endl;
    // std::cout << "omega = " << omega.transpose() << std::endl;
    // std::cout << "norm_weights_eff = " << norm_weights_eff.transpose() << std::endl;
    // return sum_weights * (omega.array().log() * norm_weights.array()).sum();
  }

  /// @param[out] omega Probability vector of length `n_obs + supp_adj`.
  /// @param[out] lambda Dual problem vector of size `n_eqs`.  
  /// @param[in] G Moment matrix of size `n_eqs x n_obs`.
  /// @param[in] weights Vector of weights of length `n_obs`.  Not assumed to be normalized.
  ///
  /// @return Value of empirical loglikelihood.
  template <class Type>
  inline Type GenEL<Type>::logel_full(RefVector_t<Type> omega,
				      RefVector_t<Type> lambda,
				      cRefMatrix_t<Type>& G,
				      cRefVector_t<Type>& weights) {
    Ref<const Matrix_t<Type> > G_eff = supp_G(G);
    Type sum_weights;
    Ref<const Vector_t<Type> > norm_weights_eff = supp_norm_weights(weights, sum_weights);
    return logel_full_impl(omega, lambda, G_eff, norm_weights_eff, sum_weights);
  }


  /// @param[in] G Moment matrix of size `n_eqs x n_obs`.
  /// @param[in] weights Vector of weights of length `n_obs`.  Not assumed to be normalized.
  /// @return Value of empirical loglikelihood.
  template <class Type>
  inline Type GenEL<Type>::logel(cRefMatrix_t<Type>& G,
				 cRefVector_t<Type>& weights) {
    Ref<const Matrix_t<Type> > G_eff = supp_G(G);
    Type sum_weights;
    Ref<const Vector_t<Type> > norm_weights_eff = supp_norm_weights(weights, sum_weights);
    // Type sum_weights = weights.sum();
    // norm_weights_.head(n_obs_) = weights;
    // if (supp_adj_) {
    //   norm_weights_(n_obs_) = weight_adj_;
    //   sum_weights += weight_adj_;
    // }
    // norm_weights_.head(n_obs2_) /= sum_weights;
    return logel_impl(G_eff, norm_weights_eff, sum_weights);
  }

  /// @param[in] G Moment matrix of size `n_eqs x n_obs`.
  /// @return The value of the empirical loglikelihood.
  template <class Type>
  inline Type GenEL<Type>::logel(cRefMatrix_t<Type>& G) {
    Ref<const Matrix_t<Type> > G_eff = supp_G(G);
    // Type sum_weights = Type(n_obs2_);
    // norm_weights_.head(n_obs2_) = Vector_t<Type> ::Ones(n_obs2_)/sum_weights;
    norm_weights_.head(n_obs_).fill(1.0);
    Type sum_weights;
    Ref<const Vector_t<Type> > norm_weights_eff = supp_norm_weights(norm_weights_.head(n_obs_), sum_weights);
    return logel_impl(G_eff, norm_weights_eff, sum_weights);
  }

  /// @param[out] dldG Matrix of size `n_eqs x n_obs` containing the gradient of logel with respect to G.
  /// @param[in] omega Probability vector of length `n_obs + supp_adj`.
  /// @param[in] lambda Dual problem vector of size `n_eqs`.
  /// @param[in] sum_weights Sum of the weights prior to normalization.
  template <class Type>
  inline void GenEL<Type>::logel_grad(RefMatrix_t<Type> dldG,
				      cRefVector_t<Type>& omega,
				      cRefVector_t<Type>& lambda,
				      Type sum_weights) {
    dldG.noalias() = lambda * omega.head(n_obs_).transpose();
    if(supp_adj_) {
      dldG.colwise() -= omega(n_obs_) * supp_a_/n_obs_ * lambda;
    }
    dldG *= sum_weights;
    return;
  }

  /// @param[out] dldG Matrix of size `n_eqs x n_obs` containing the gradient of logel with respect to G.
  /// @param[in] G Moment matrix of size `n_eqs x n_obs`.
  /// @return For convenience, returns the value of the empirical loglikelihood.
  template <class Type>
  inline Type GenEL<Type>::logel_grad(RefMatrix_t<Type> dldG,
				      cRefMatrix_t<Type>& G,
				      cRefVector_t<Type>& weights) {
    // Type sum_weights = weights.sum();
    // if (supp_adj_) {
    //   sum_weights += weight_adj_;
    // }
    // Type ll = logel(G, weights);
    Ref<const Matrix_t<Type> > G_eff = supp_G(G);
    Type sum_weights;
    Ref<const Vector_t<Type> > norm_weights_eff = supp_norm_weights(weights, sum_weights);
    Type ll = logel_impl(G_eff, norm_weights_eff, sum_weights);    
    // use internal values of omega_ and lambda_
    logel_grad(dldG, omega_.head(n_obs2_), lambda_, sum_weights);
    return ll;
  }

  /// @param[out] dldG Matrix of size `n_eqs x n_obs` containing the gradient of logel with respect to G.
  /// @param[in] G Moment matrix of size `n_eqs x n_obs`.
  /// @return For convenience, returns the value of the empirical loglikelihood.
  template <class Type>
  inline Type GenEL<Type>::logel_grad(RefMatrix_t<Type> dldG,
				      cRefMatrix_t<Type>& G) {
    Ref<const Matrix_t<Type> > G_eff = supp_G(G);
    norm_weights_.head(n_obs_).fill(1.0);
    Type sum_weights;
    Ref<const Vector_t<Type> > norm_weights_eff = supp_norm_weights(norm_weights_.head(n_obs_), sum_weights);
    Type ll = logel_impl(G_eff, norm_weights_eff, sum_weights);
    // Type ll = logel(G);
    // use internal values of omega_ and lambda_
    logel_grad(dldG, omega_.head(n_obs2_), lambda_, sum_weights);
    return ll;
  }
    
} // end namespace flexEL

#endif
