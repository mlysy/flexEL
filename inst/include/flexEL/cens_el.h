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

  class CensEL {
    
  private:
    
    // dim of G matrix
    int n_obs_; // number of columns
    int n_eqs_; // number of rows    
    int n_obs1_; // n_obs1_ = n_obs_+1: for memory allocation
    bool supp_adj_; // support adjustment flag
    double smooth_s_; // continuity correction parameter
    
    // internal variables for EM algorithm
    VectorXd lambda_; // NR dual vector
    VectorXd weights_; // EM weights
    VectorXd omega_; // probability vector
    VectorXd omega_tilde_; // for E-step
    VectorXd delta_aug_; // censoring indicators
    VectorXd epsilon_aug_; // residuals
    int max_iter_; // maximum number of iterations
    int em_iter_; // actual number of iterations
    double abs_tol_; // absolute tolerance
    double em_err_; // actual absolute error
    bool has_converged_nr_; // flag for NR convergence failure

    // implementation functions
    void omega_hat_impl(Ref<VectorXd> omega,
			const Ref<const MatrixXd>& G,
			const Ref<const VectorXd>& delta,
			const Ref<const VectorXd>& epsilon,
			GenEL& gel);
    double logel_omega_impl(const Ref<const VectorXd>& omega,
			    const Ref<const VectorXd>& delta,
			    const Ref<const VectorXd>& epsilon);
    
  public:
    
    /// Constructor.
    CensEL(int n_obs, int n_eqs);
    /// Set the maximum number of EM iterations.
    void set_max_iter(int max_iter);
    /// Set the absolute tolerance for the EM algorithm.
    void set_abs_tol(double abs_tol);
    /// Set the smoothing parameter.
    void set_smooth(double s);
    /// Get number of observations.
    int get_n_obs();
    /// Get number of estimating equations.
    int get_n_eqs();
    /// Get the diagnostics for the EM run.
    void get_diag(int& em_iter, double& em_err);
    /// Get the maximum number of EM iterations.
    int get_max_iter();
    /// Get the absolute tolerance of the EM algorithm.
    double get_abs_tol();
    /// Get the smoothing parameter.
    double get_smooth();

    /// Check whether EM algorithm has converged.
    bool has_converged_em();

    /// Censoring vector with support adjustment.
    Ref<const VectorXd> supp_delta(const Ref<const VectorXd>& delta);
    /// Residual vector with support adjustment.
    Ref<const VectorXd> supp_epsilon(const Ref<const VectorXd>& epsilon);
    
    /// Calculate the expected weights in the EM algorithm.
    void expected_weights(Ref<VectorXd> weights,
			  const Ref<const VectorXd>& delta,
			  const Ref<const VectorXd>& epsilon,
			  const Ref<const VectorXd>& omega);
        
    /// Calculate the profile probability weights.
    void omega_hat(Ref<VectorXd> omega,
                   const Ref<const MatrixXd>& G,
                   const Ref<const VectorXd>& delta,
                   const Ref<const VectorXd>& epsilon,
		   GenEL& gel);

    /// Calculate the censored empirical loglikelihood given the probability weights.
    double logel_omega(const Ref<const VectorXd>& omega,
		       const Ref<const VectorXd>& delta,
		       const Ref<const VectorXd>& epsilon);
    /// Calculate the censored empirical loglikelihood.
    double logel(const Ref<const MatrixXd>& G,
                 const Ref<const VectorXd>& delta,
                 const Ref<const VectorXd>& epsilon,
		 GenEL& gel);
    
  };
  
  /// @param[in] n_obs    Number of observations.
  /// @param[in] n_eqs    Number of estimating equations.
  inline CensEL::CensEL(int n_obs, int n_eqs) {
    // assign internal values
    n_obs_ = n_obs;
    n_eqs_ = n_eqs;
    n_obs1_ = n_obs_+1; // space for augmented delta and epsilon
    supp_adj_ = false;
    has_converged_nr_ = false;
    // memory allocation
    lambda_ = VectorXd::Zero(n_eqs_);
    weights_ = VectorXd::Zero(n_obs1_);
    delta_aug_ = VectorXd::Zero(n_obs1_);
    epsilon_aug_ = VectorXd::Zero(n_obs1_);
    omega_ = VectorXd::Zero(n_obs1_);
    omega_tilde_ = VectorXd::Zero(n_obs1_);
    // support adjustment dummy values
    delta_aug_(n_obs_) = 0;
    epsilon_aug_(n_obs_) = -std::numeric_limits<double>::infinity();
    // default tuning parameters
    em_iter_ = 0;
    em_err_ = 0.0;
    set_smooth(10.0);
    set_max_iter(100);
    set_abs_tol(1e-3);
  }
  
  /// @param[in] max_iter Maximum number of EM iterations.
  inline void CensEL::set_max_iter(int max_iter) {
    max_iter_ = max_iter;
    return;
  }
  
  /// @param[in] abs_tol Relative tolerance for the EM algorithm.
  inline void CensEL::set_abs_tol(double abs_tol) {
    abs_tol_ = abs_tol;
    return;
  }
  
  /// @param[in] s      Tuning parameter for continuity correction.
  inline void CensEL::set_smooth(double s) {
    smooth_s_ = s; 
    return;
  }

  inline int CensEL::get_n_obs() {
    return n_obs_;
  }
  
  inline int CensEL::get_n_eqs() {
    return n_eqs_;
  }

  
  /// @param[out] em_iter Number of EM iterations.
  /// @param[out] em_err Maximum relative difference between the value of logel in the last two EM steps.
  inline void CensEL::get_diag(int& em_iter, double& em_err) {
    em_iter = em_iter_;
    em_err = em_err_;
    return;
  }

  inline double CensEL::get_smooth() {
    return smooth_s_;
  }
  
  inline int CensEL::get_max_iter() {
    return max_iter_;
  }
  
  inline double CensEL::get_abs_tol() {
    return abs_tol_;
  }
  
  /// @return `true` if (i) all the Newton-Raphson M-steps converged and (ii) the desired error tolerance `abs_tol` has been reached in less than `max_iter` EM steps.
  inline bool CensEL::has_converged_em() {
    return has_converged_nr_ &&
      (!((em_err_ > abs_tol_) && (em_iter_ == max_iter_)));
  }


  /// @param[in] delta Censoring vector of length `n_obs` or `n_obs+1`. If `supp_adj = false`, the former is required.  If `supp_adj = true` and the former is provided, support adjustment is performed.  If `supp_adj = true` and `delta.size() == n_obs+1`, assumes that support has already been adjusted.
  /// return A reference to `delta` if no support adjustment is necessary (not needed or already performed); otherwise a reference to `delta_aug_`. 
  Ref<const VectorXd> CensEL::supp_delta(const Ref<const VectorXd>& delta) {
    if(supp_adj_ && delta.size() == n_obs_) {
      delta_aug_.head(n_obs_) = delta;
      return Ref<const VectorXd>(delta_aug_);
    } else {
      return Ref<const VectorXd>(delta);
    }
  }
  
  /// @param[in] epsilon Residual vector of length `n_obs` or `n_obs+1`. If `supp_adj = false`, the former is required.  If `supp_adj = true` and the former is provided, support adjustment is performed.  If `supp_adj = true` and `epsilon.size() == n_obs+1`, assumes that support has already been adjusted.
  /// return A reference to `epsilon` if no support adjustment is necessary (not needed or already performed); otherwise a reference to `epsilon_aug_`. 
  Ref<const VectorXd> CensEL::supp_epsilon(const Ref<const VectorXd>& epsilon) {
    if(supp_adj_ && epsilon.size() == n_obs_) {
      epsilon_aug_.head(n_obs_) = epsilon;
      return Ref<const VectorXd>(epsilon_aug_);
    } else {
      return Ref<const VectorXd>(epsilon);
    }
  }


  /// @param[out] weights Weight vector of length `n_obs` or `n_obs+1`.
  /// @param[in] delta Censoring vector of the same length as `weights`.
  /// @param[in] epsilon Residual vector of the same length as `weights`.
  /// @param[in] omega Probability vector of the same length as `weights`.
  inline void CensEL::expected_weights(Ref<VectorXd> weights,
				       const Ref<const VectorXd>& delta,
				       const Ref<const VectorXd>& epsilon,
				       const Ref<const VectorXd>& omega) {
    int n_obs2 = weights.size();
    weights.setZero();
    double denom;
    for(int ii=0; ii<n_obs2; ii++) {
      if(delta(ii) > 0.5) {
        // avoid checking double > 0.0: delta should only ever be 0 or 1!
	weights(ii) += 1.0;
      } else {
	smooth_indicator(omega_tilde_.head(n_obs2), epsilon(ii), epsilon, smooth_s_);
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
  /// @param[in] gel Reference to `GenEL` object which performs the M-step of the EM algorithm.
  inline void CensEL::omega_hat_impl(Ref<VectorXd> omega,
				     const Ref<const MatrixXd>& G,
				     const Ref<const VectorXd>& delta,
				     const Ref<const VectorXd>& epsilon,
				     GenEL& gel) {
    int n_obs2 = omega.size();
    // initialize quantities
    omega.fill(1.0/double(n_obs2));
    expected_weights(weights_.head(n_obs2), delta, epsilon, omega);
    double logel_old = gel.logel(G);
    double logel;
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
      if(em_err_ < abs_tol_) break;
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
  inline double CensEL::logel_omega_impl(const Ref<const VectorXd>& omega,
					 const Ref<const VectorXd>& delta,
					 const Ref<const VectorXd>& epsilon) {
    int n_obs2 = omega.size();
    double logel = 0.0;
    for(int ii=0; ii<n_obs2; ii++) {
      if(delta(ii)) {
	logel += log(omega(ii));
      } else {
	smooth_indicator(omega_tilde_.head(n_obs2), epsilon(ii), epsilon, smooth_s_);
	logel += log((omega_tilde_.head(n_obs2).array() * omega.array()).sum()); 
      }
    }
    return logel;
  }

  /// @param[out] omega  Probability vector of length `n_obs + supp_adj`.
  /// @param[in] G       Moment matrix of size `n_eqs x n_obs` or `n_eqs x (n_obs + supp_adj)`.  If `supp_adj = false`, the former is required.  If `supp_adj = true` and the former is provided, support adjustment is performed.  If `supp_adj = true` and `G.cols() == n_obs + 1`, assumes that support has already been corrected. 
  /// @param[in] delta   Vector of censoring indicators of length `n_obs` or `n_obs + supp_adj` following same logic as above.
  /// @param[in] gel Reference to `GenEL` object which performs the M-step of the EM algorithm.
  inline void CensEL::omega_hat(Ref<VectorXd> omega,
                                const Ref<const MatrixXd>& G,
                                const Ref<const VectorXd>& delta,
                                const Ref<const VectorXd>& epsilon,
				GenEL& gel) {
    supp_adj_ = gel.get_supp_adj();
    Ref<const VectorXd> delta_eff = supp_delta(delta);
    Ref<const VectorXd> epsilon_eff = supp_epsilon(epsilon);
    Ref<const MatrixXd> G_eff = gel.supp_G(G);
    omega_hat_impl(omega, G_eff, delta_eff, epsilon_eff, gel);
    return;
  }

  /// @param[in] omega Probability vector of length `n_obs + supp_adj`.
  /// @param[in] delta Censoring vector of length `n_obs` or `n_obs + `supp_adj`.
  /// @param[in] epsilon Residual vector of length `n_obs` or `n_obs + supp_adj`.
  /// @return Value of the censored empirical loglikelihood.
  inline double CensEL::logel_omega(const Ref<const VectorXd>& omega,
				    const Ref<const VectorXd>& delta,
				    const Ref<const VectorXd>& epsilon) {
    supp_adj_ = omega.size() == n_obs1_;
    Ref<const VectorXd> delta_eff = supp_delta(delta);
    Ref<const VectorXd> epsilon_eff = supp_epsilon(epsilon);
    return logel_omega_impl(omega, delta_eff, epsilon_eff);
  }
  
  /// @param[in] G       Moment matrix of size `n_eqs x n_obs`. 
  /// @param[in] delta   Censoring indicator vector of length `n_obs`.
  /// @param[in] epsilon Residual vector of length `n_obs`.
  /// @param[in] gel     Reference to `GenEL` object which performs the M-step of the EM algorithm.
  /// @return Value of the censored empirical loglikelihood.
  inline double CensEL::logel(const Ref<const MatrixXd>& G,
                              const Ref<const VectorXd>& delta,
                              const Ref<const VectorXd>& epsilon,
			      GenEL& gel) {
    supp_adj_ = gel.get_supp_adj();
    int n_obs2 = n_obs_ + supp_adj_;
    Ref<const VectorXd> delta_eff = supp_delta(delta);
    Ref<const VectorXd> epsilon_eff = supp_epsilon(epsilon);
    Ref<const MatrixXd> G_eff = gel.supp_G(G);
    omega_hat_impl(omega_.head(n_obs2), G_eff, delta_eff, epsilon_eff, gel);
    return logel_omega_impl(omega_.head(n_obs2), delta_eff, epsilon_eff);
  }
  
}

#endif
