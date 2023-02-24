/// @file cens_el.h

#ifndef FLEXEL_CENS_EL_H
#define FLEXEL_CENS_EL_H

#include <Eigen/Dense>
#include "utils.h"
#include "gen_el.h"
#include <limits> // for std::numeric_limits
#include <iostream>
// #include "sort_order.h"
// #include "ind_smooth.h" 

namespace flexEL {

  using namespace Eigen;

  template <class Type>
  class CensEL {
    
  private:
    
    // dim of G matrix
    int n_obs_; // number of columns
    int n_eqs_; // number of rows    
    int n_obs1_; // n_obs1_ = n_obs_+1: for memory allocation
    bool supp_adj_; // support adjustment flag
    Type smooth_s_; // continuity correction parameter
    
    // internal variables for EM algorithm
    Vector_t<Type> lambda_; // NR dual vector
    Vector_t<Type> weights_; // EM weights
    Vector_t<Type> omega_; // probability vector
    Vector_t<Type> omega_tilde_; // for E-step
    Vector_t<Type> delta_aug_; // censoring indicators
    Vector_t<Type> epsilon_aug_; // residuals
    int max_iter_; // maximum number of iterations
    int em_iter_; // actual number of iterations
    Type abs_tol_; // absolute tolerance
    Type em_err_; // actual absolute error
    bool has_converged_nr_; // flag for NR convergence failure

    // implementation functions
    void omega_hat_impl(RefVector_t<Type> omega,
			cRefMatrix_t<Type>& G,
			cRefVector_t<Type>& delta,
			cRefVector_t<Type>& epsilon,
			GenEL<Type>& gel);
    Type logel_omega_impl(cRefVector_t<Type>& omega,
			  cRefVector_t<Type>& delta,
			  cRefVector_t<Type>& epsilon);
    
  public:
    
    /// Constructor.
    CensEL(int n_obs, int n_eqs);
    /// Set the maximum number of EM iterations.
    void set_max_iter(int max_iter);
    /// Set the absolute tolerance for the EM algorithm.
    void set_abs_tol(Type abs_tol);
    /// Set the smoothing parameter.
    void set_smooth(Type s);
    /// Get number of observations.
    int get_n_obs();
    /// Get number of estimating equations.
    int get_n_eqs();
    /// Get the diagnostics for the EM run.
    void get_diag(int& em_iter, Type& em_err);
    /// Get the maximum number of EM iterations.
    int get_max_iter();
    /// Get the absolute tolerance of the EM algorithm.
    Type get_abs_tol();
    /// Get the smoothing parameter.
    Type get_smooth();

    /// Check whether EM algorithm has converged.
    bool has_converged_em();

    /// Censoring vector with support adjustment.
    Ref<const Vector_t<Type> > supp_delta(cRefVector_t<Type>& delta);
    /// Residual vector with support adjustment.
    Ref<const Vector_t<Type> > supp_epsilon(cRefVector_t<Type>& epsilon);
    
    /// Calculate the expected weights in the EM algorithm.
    void expected_weights(RefVector_t<Type> weights,
			  cRefVector_t<Type>& delta,
			  cRefVector_t<Type>& epsilon,
			  cRefVector_t<Type>& omega);
        
    /// Calculate the profile probability weights.
    void omega_hat(RefVector_t<Type> omega,
                   cRefMatrix_t<Type>& G,
                   cRefVector_t<Type>& delta,
                   cRefVector_t<Type>& epsilon,
		   GenEL<Type>& gel);

    /// Calculate the censored empirical loglikelihood given the probability weights.
    Type logel_omega(cRefVector_t<Type>& omega,
		     cRefVector_t<Type>& delta,
		     cRefVector_t<Type>& epsilon);
    /// Calculate the censored empirical loglikelihood.
    Type logel(cRefMatrix_t<Type>& G,
	       cRefVector_t<Type>& delta,
	       cRefVector_t<Type>& epsilon,
	       GenEL<Type>& gel);
    
  };
  
  /// @param[in] n_obs    Number of observations.
  /// @param[in] n_eqs    Number of estimating equations.
  template <class Type>
  inline CensEL<Type>::CensEL(int n_obs, int n_eqs) {
    // assign internal values
    n_obs_ = n_obs;
    n_eqs_ = n_eqs;
    n_obs1_ = n_obs_+1; // space for augmented delta and epsilon
    supp_adj_ = false;
    has_converged_nr_ = false;
    // memory allocation
    lambda_ = Vector_t<Type>::Zero(n_eqs_);
    weights_ = Vector_t<Type>::Zero(n_obs1_);
    delta_aug_ = Vector_t<Type>::Zero(n_obs1_);
    epsilon_aug_ = Vector_t<Type>::Zero(n_obs1_);
    omega_ = Vector_t<Type>::Zero(n_obs1_);
    omega_tilde_ = Vector_t<Type>::Zero(n_obs1_);
    // support adjustment dummy values
    delta_aug_(n_obs_) = 0;
    epsilon_aug_(n_obs_) = -std::numeric_limits<Type>::infinity();
    // default tuning parameters
    em_iter_ = 0;
    em_err_ = 0.0;
    set_smooth(10.0);
    set_max_iter(100);
    set_abs_tol(1e-3);
  }
  
  /// @param[in] max_iter Maximum number of EM iterations.
  template <class Type>
  inline void CensEL<Type>::set_max_iter(int max_iter) {
    max_iter_ = max_iter;
    return;
  }
  
  /// @param[in] abs_tol Relative tolerance for the EM algorithm.
  template <class Type>
  inline void CensEL<Type>::set_abs_tol(Type abs_tol) {
    abs_tol_ = abs_tol;
    return;
  }
  
  /// @param[in] s      Tuning parameter for continuity correction.
  template <class Type>
  inline void CensEL<Type>::set_smooth(Type s) {
    smooth_s_ = s; 
    return;
  }

  template <class Type>
  inline int CensEL<Type>::get_n_obs() {
    return n_obs_;
  }
  
  template <class Type>
  inline int CensEL<Type>::get_n_eqs() {
    return n_eqs_;
  }

  
  /// @param[out] em_iter Number of EM iterations.
  /// @param[out] em_err Maximum relative difference between the value of logel in the last two EM steps.
  template <class Type>
  inline void CensEL<Type>::get_diag(int& em_iter, Type& em_err) {
    em_iter = em_iter_;
    em_err = em_err_;
    return;
  }

  template <class Type>
  inline Type CensEL<Type>::get_smooth() {
    return smooth_s_;
  }
  
  template <class Type>
  inline int CensEL<Type>::get_max_iter() {
    return max_iter_;
  }
  
  template <class Type>
  inline Type CensEL<Type>::get_abs_tol() {
    return abs_tol_;
  }
  
  /// @return `true` if (i) all the Newton-Raphson M-steps converged and (ii) the desired error tolerance `abs_tol` has been reached in less than `max_iter` EM steps.
  template <class Type>
  inline bool CensEL<Type>::has_converged_em() {
    return has_converged_nr_ &&
      (!((em_err_ > abs_tol_) && (em_iter_ == max_iter_)));
  }


  /// @param[in] delta Censoring vector of length `n_obs` or `n_obs+1`. If `supp_adj = false`, the former is required.  If `supp_adj = true` and the former is provided, support adjustment is performed.  If `supp_adj = true` and `delta.size() == n_obs+1`, assumes that support has already been adjusted.
  /// return A reference to `delta` if no support adjustment is necessary (not needed or already performed); otherwise a reference to `delta_aug_`. 
  template <class Type>
  Ref<const Vector_t<Type> > CensEL<Type>::supp_delta(cRefVector_t<Type>& delta) {
    if(supp_adj_ && delta.size() == n_obs_) {
      delta_aug_.head(n_obs_) = delta;
      return Ref<const Vector_t<Type> >(delta_aug_);
    } else {
      return Ref<const Vector_t<Type> >(delta);
    }
  }
  
  /// @param[in] epsilon Residual vector of length `n_obs` or `n_obs+1`. If `supp_adj = false`, the former is required.  If `supp_adj = true` and the former is provided, support adjustment is performed.  If `supp_adj = true` and `epsilon.size() == n_obs+1`, assumes that support has already been adjusted.
  /// return A reference to `epsilon` if no support adjustment is necessary (not needed or already performed); otherwise a reference to `epsilon_aug_`. 
  template <class Type>
  Ref<const Vector_t<Type> > CensEL<Type>::supp_epsilon(cRefVector_t<Type>& epsilon) {
    if(supp_adj_ && epsilon.size() == n_obs_) {
      epsilon_aug_.head(n_obs_) = epsilon;
      return Ref<const Vector_t<Type> >(epsilon_aug_);
    } else {
      return Ref<const Vector_t<Type> >(epsilon);
    }
  }


  /// @param[out] weights Weight vector of length `n_obs` or `n_obs+1`.
  /// @param[in] delta Censoring vector of the same length as `weights`.
  /// @param[in] epsilon Residual vector of the same length as `weights`.
  /// @param[in] omega Probability vector of the same length as `weights`.
  template <class Type>
  inline void CensEL<Type>::expected_weights(RefVector_t<Type> weights,
					     cRefVector_t<Type>& delta,
					     cRefVector_t<Type>& epsilon,
					     cRefVector_t<Type>& omega) {
    int n_obs2 = weights.size();
    weights.setZero();
    Type denom;
    for(int ii=0; ii<n_obs2; ii++) {
      if(delta(ii) > 0.5) {
        // avoid checking Type > 0.0: delta should only ever be 0 or 1!
	weights(ii) += 1.0;
      } else {
	smooth_indicator<Type>(omega_tilde_.head(n_obs2), epsilon(ii), epsilon, smooth_s_);
	if(ii == n_obs_) {
	  // convert 1/(1 + exp(Inf-Inf)) to 0.5
	  omega_tilde_(ii) = 0.5;
	}
	omega_tilde_.head(n_obs2).array() *= omega.array();
	denom = (omega_tilde_.head(n_obs2).array()).sum();
	weights.array() += omega_tilde_.head(n_obs2).array() / denom;
	// omega_tilde_.head(n_obs2).array() /= denom;
	// weights(ii) = delta(ii) + ((1.0 - delta.array()) * omega_tilde_.head(n_obs2).array()).sum();
      }
    }
    return;
  }

  /// Convergence of the algorithm is determined as follows:
  ///
  /// - After each M-step, convergence of the NR algorithm is checked.  If it failed, the EM algorithm is deemed not to have converged.  However, it does not stop running.
  ///
  /// - The EM algorithm stops when the termination criterion is reached, i.e., either `max_iter` has been reached or the difference between `logel` values is less than `abs_tol`.
  ///
  /// - The algorithm can also stop when `gel.logel()` returns a `nan`.
  /// Once the EM algorithm stops, we check whether every NR algorithm has converged and also whether the absolute tolerance has been reached.  If these are both true, the EM algorithm is deemed to have converged.
  ///
  /// @param[out] omega  Probability vector of length `n_obs` or `n_obs+1`.
  /// @param[in] G       Moment matrix with `n_eqs` rows and `omega.size()` columns. 
  /// @param[in] delta   Vector of censoring indicators of length `omega.size()`.
  /// @param[in] gel Reference to `GenEL<Type>` object which performs the M-step of the EM algorithm.
  template <class Type>
  inline void CensEL<Type>::omega_hat_impl(RefVector_t<Type> omega,
					   cRefMatrix_t<Type>& G,
					   cRefVector_t<Type>& delta,
					   cRefVector_t<Type>& epsilon,
					   GenEL<Type>& gel) {
    int n_obs2 = omega.size();
    // initialize quantities
    omega.fill(1.0/Type(n_obs2));
    expected_weights(weights_.head(n_obs2), delta, epsilon, omega);
    Type logel_old = gel.logel(G);
    Type logel;
    int ii;
    has_converged_nr_ = true;
    // printf("max_iter_nr = %i, rel_tol_nr = %f\n",
    // 	   gel.get_max_iter(), gel.get_rel_tol());
    // std::cout << "weights:" << std::endl << weights_.head(n_obs2).transpose() << std::endl;
    for(ii=0; ii<max_iter_; ii++) {
      // std::cout << "step " << ii << std::endl;
      // M-step
      gel.lambda_nr(lambda_, G, weights_.head(n_obs2));
      gel.omega_hat(omega, lambda_, G, weights_.head(n_obs2));
      logel = gel.logel_omega(omega, weights_.head(n_obs2));
      // std::cout << "lambda:" << std::endl << lambda_.transpose() << std::endl;
      // std::cout << "omega:" << std::endl << omega.transpose() << std::endl;
      // printf("logel = %f, logel_old = %f\n", logel, logel_old);
      em_err_ = abs(logel - logel_old);
      // printf("em_err = %f, abs_tol = %f\n", em_err_, abs_tol_);
      // printf("gel.has_converged_nr() = %i\n", gel.has_converged_nr());
      if(has_converged_nr_ && !gel.has_converged_nr()) {
	has_converged_nr_ = false;
      }
      if(isnan(em_err_) || em_err_ < abs_tol_) break;
      logel_old = logel;
      gel.set_lambda0(lambda_); // update starting point
      // E-step
      expected_weights(weights_.head(n_obs2), delta, epsilon, omega);
      // std::cout << "weights:" << std::endl << weights_.head(n_obs2).transpose() << std::endl;
    }
    em_iter_ = ii;
    return;
  }

  /// @param[in] omega Probability vector of length `n_obs + supp_adj`.
  /// @param[in] delta Censoring vector of same length as `omega`.
  /// @param[in] epsilon Residual vector of same length as `omega`.
  /// @return Value of the censored empirical loglikelihood.
  template <class Type>
  inline Type CensEL<Type>::logel_omega_impl(cRefVector_t<Type>& omega,
					     cRefVector_t<Type>& delta,
					     cRefVector_t<Type>& epsilon) {
    int n_obs2 = omega.size();
    Type logel = 0.0;
    for(int ii=0; ii<n_obs2; ii++) {
      if(delta(ii)) {
	logel += log(omega(ii));
      } else {
	smooth_indicator<Type>(omega_tilde_.head(n_obs2), epsilon(ii), epsilon, smooth_s_);
	logel += log((omega_tilde_.head(n_obs2).array() * omega.array()).sum()); 
      }
    }
    return logel;
  }

  /// @param[out] omega  Probability vector of length `n_obs + supp_adj`.
  /// @param[in] G       Moment matrix of size `n_eqs x n_obs` or `n_eqs x (n_obs + supp_adj)`.  If `supp_adj = false`, the former is required.  If `supp_adj = true` and the former is provided, support adjustment is performed.  If `supp_adj = true` and `G.cols() == n_obs + 1`, assumes that support has already been corrected. 
  /// @param[in] delta   Vector of censoring indicators of length `n_obs` or `n_obs + supp_adj` following same logic as above.
  /// @param[in] gel Reference to `GenEL<Type>` object which performs the M-step of the EM algorithm.
  template <class Type>
  inline void CensEL<Type>::omega_hat(RefVector_t<Type> omega,
				      cRefMatrix_t<Type>& G,
				      cRefVector_t<Type>& delta,
				      cRefVector_t<Type>& epsilon,
				      GenEL<Type>& gel) {
    supp_adj_ = gel.get_supp_adj();
    Ref<const Vector_t<Type> > delta_eff = supp_delta(delta);
    Ref<const Vector_t<Type> > epsilon_eff = supp_epsilon(epsilon);
    Ref<const Matrix_t<Type> > G_eff = gel.supp_G(G);
    omega_hat_impl(omega, G_eff, delta_eff, epsilon_eff, gel);
    return;
  }

  /// @param[in] omega Probability vector of length `n_obs + supp_adj`.
  /// @param[in] delta Censoring vector of length `n_obs` or `n_obs + `supp_adj`.
  /// @param[in] epsilon Residual vector of length `n_obs` or `n_obs + supp_adj`.
  /// @return Value of the censored empirical loglikelihood.
  template <class Type>
  inline Type CensEL<Type>::logel_omega(cRefVector_t<Type>& omega,
					cRefVector_t<Type>& delta,
					cRefVector_t<Type>& epsilon) {
    supp_adj_ = omega.size() == n_obs1_;
    Ref<const Vector_t<Type> > delta_eff = supp_delta(delta);
    Ref<const Vector_t<Type> > epsilon_eff = supp_epsilon(epsilon);
    return logel_omega_impl(omega, delta_eff, epsilon_eff);
  }
  
  /// @param[in] G       Moment matrix of size `n_eqs x n_obs`. 
  /// @param[in] delta   Censoring indicator vector of length `n_obs`.
  /// @param[in] epsilon Residual vector of length `n_obs`.
  /// @param[in] gel     Reference to `GenEL<Type>` object which performs the M-step of the EM algorithm.
  /// @return Value of the censored empirical loglikelihood.
  template <class Type>
  inline Type CensEL<Type>::logel(cRefMatrix_t<Type>& G,
				  cRefVector_t<Type>& delta,
				  cRefVector_t<Type>& epsilon,
				  GenEL<Type>& gel) {
    supp_adj_ = gel.get_supp_adj();
    int n_obs2 = n_obs_ + supp_adj_;
    Ref<const Vector_t<Type> > delta_eff = supp_delta(delta);
    Ref<const Vector_t<Type> > epsilon_eff = supp_epsilon(epsilon);
    Ref<const Matrix_t<Type> > G_eff = gel.supp_G(G);
    omega_hat_impl(omega_.head(n_obs2), G_eff, delta_eff, epsilon_eff, gel);
    return logel_omega_impl(omega_.head(n_obs2), delta_eff, epsilon_eff);
  }
  
}

#endif
