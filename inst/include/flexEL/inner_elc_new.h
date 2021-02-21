/// @file inner_elc.h

#ifndef INNER_ELC_H
#define INNER_ELC_H

#include <RcppEigen.h>
#include "inner_el.h"
#include "sort_order.h"

namespace flexEL {

  using namespace Eigen;

  class CensEL {
    
  private:
    
    GenEL GEL;
    
    // dimensions of G
    int n_obs_; // number of columns
    int n_eqs_; // number of rows    
    // for support modification
    bool supp_adj_;
    double supp_a_;
    int n_obs1_; // n_obs1_ = n_obs_+1: for memory allocation
    int n_obs2_; // n_obs2_ = n_obs_+support: used n_obs_ in calculation
    
    // internal variables for EM algorithm
    VectorXd omega_init_; // initial value for omegas, in case of reset
    VectorXd delta_; // censoring indicators
    VectorXd epsilon_;
    VectorXi eps_ord_; 
    VectorXd norm_weights_;
    VectorXd psot_; // partial sum of omegatildas

    // internals for EM options
    int max_iter_em_;
    double abs_tol_; // absolute tolerance
    
    // helper function for eval_weights: calculate partial sum of omegas
    void eval_pso(double &psos,
                  const int ii,
                  const Ref<const VectorXd>& omega); // Partial sum of omega_jj s.t. eps_jj >= eps_ii
    
    
  public:
    
    /// Constructor.
    CensEL(int n_obs, int n_eqs);
    /// Set the maximum number of Newton-Raphson iterations.
    void set_max_iter_nr(int max_iter);
    /// Set the maximum number of EM iterations.
    void set_max_iter_em(int max_iter);
    /// Set the relative tolerance for the Newton-Raphson algorithm.
    void set_rel_tol(double rel_tol);
    /// Set the absolute tolerance for the EM algorithm.
    void set_abs_tol(double abs_tol);
    /// Set the support adjustment flag.
    void set_supp_adj(bool supp_adj, double a);
    void set_supp_adj(bool supp_adj);
    
    /// Calculate weights according to epsilons
    void eval_weights(Ref<VectorXd> weights,
                      const Ref<const VectorXd>& delta,
                      const Ref<const VectorXd>& epsilon,
                      const Ref<const VectorXd>& omega);
    
    void omega_hat(Ref<VectorXd> omega,
                   const Ref<const MatrixXd>& G,
                   const Ref<const VectorXd>& delta,
                   const Ref<const VectorXd>& epsilon);
    
    double logel_omega(const Ref<const VectorXd>& delta,
                       const Ref<const VectorXd>& epsilon,
                       const Ref<const VectorXd>& omega);
    
  };
  
  /// @param[in] n_obs    Number of observations.
  /// @param[in] n_eqs    Number of estimating equations.
  inline CensEL::CensEL(int n_obs, int n_eqs) : GEL(n_obs, n_eqs) {
    
    // assign internal values
    n_obs_ = n_obs;
    n_eqs_ = n_eqs;
    supp_adj_ = false; // no support correction by default
    n_obs1_ = n_obs_+1; // space for augmented G
    n_obs2_ = n_obs_ + supp_adj_;

    omega_init_ = VectorXd::Zero(n_obs1_).array() + 1.0/(double)n_obs1_; 
    delta_ = VectorXd::Zero(n_obs1_);
    epsilon_ = VectorXd::Zero(n_obs1_);
    eps_ord_ = VectorXi::Zero(n_obs1_);
    norm_weights_ = VectorXd::Zero(n_obs1_);
    psot_ = VectorXd::Zero(n_obs1_);
  }
  
  /// @param[in] max_iter Maximum number of Newton-Raphson iterations.
  inline void CensEL::set_max_iter_nr(int max_iter) {
    GEL.set_max_iter(max_iter);
    return;
  }
  
  /// @param[in] rel_tol Relative tolerance for the Newton-Raphson algorithm.
  inline void CensEL::set_rel_tol(double rel_tol) {
    GEL.set_rel_tol(rel_tol);
    return;
  }
  
  /// @param[in] max_iter Maximum number of EM iterations.
  inline void CensEL::set_max_iter_em(int max_iter) {
    max_iter_em_ = max_iter;
    return;
  }
  
  /// @param[in] rel_tol Relative tolerance for the EM algorithm.
  inline void CensEL::set_abs_tol(double abs_tol) {
    abs_tol_ = abs_tol;
    return;
  }
  
  inline void CensEL::set_supp_adj(bool supp_adj, double a) {
    supp_adj_ = supp_adj; 
    supp_a_ = a; 
    GEL.set_supp_adj(supp_adj_, supp_a_);
    if (supp_adj) {
      epsilon_.tail(1)(0) = -INFINITY;
      delta_.tail(1)(0) = 0;
    }
    return;
  }
  
  inline void CensEL::set_supp_adj(bool supp_adj) {
    set_supp_adj(supp_adj, std::max(1.0,0.5*log(n_obs_)));
    return;
  }
  
  /// @param[out] psos Partial sum of omega.
  /// @param[in] ii Index for ii-th largest epsilon.
  /// @param[in] omega Probability vector of length `n_obs + supp_adj`.
  inline void CensEL::eval_pso(double &psos,
                                 const int ii,
                                 const Ref<const VectorXd>& omega) {
    int kk;
    for (int jj=n_obs2_-1; jj>=0; jj--) {
      // Note: eps_ord_ corresponds to epsilons in ascending order
      kk = eps_ord_(jj); // kk is the index of the jj-th largest epsilon
      psos += omega(kk);
      if (kk == ii) break; // until (backwardly) reaching ii-th largest epsilon 
    }
    return;
  }
  
  inline void CensEL::eval_weights(Ref<VectorXd> weights,
                                   const Ref<const VectorXd>& delta,
                                   const Ref<const VectorXd>& epsilon,
                                   const Ref<const VectorXd>& omega) {
    // find the indices for increasing order of epsilons 
    psot_.fill(0.0);
    int kk;
    double psos;
    delta_.head(n_obs_) = delta; 
    epsilon_.head(n_obs_) = epsilon;
    VectorXi eps_ord_ = sort_inds(epsilon_); 
    for (int ii=0; ii<n_obs2_; ii++) {
      for (int jj=0; jj<n_obs2_; jj++) {
        kk = eps_ord_(jj);
        if (delta_(kk) == 0) {
          eval_pso(psos, kk, omega);
          // to prevent dividing by 0
          if (abs(psos) >= 1e-10) psot_(ii) += omega(ii)/psos;
          // else if (omega_(ii) >= 1e-10 && EvalPSO(kk) < 1e-10) {
          //   // TODO: this means a problem
          //   // std::cout << "EvalWeights: dividing by 0 problem." << std::endl;
          // }
        }
        if (kk == ii) break;
      }
    }
    weights.array() = delta_.array() + psot_.array();
  }
  
  inline void CensEL::omega_hat(Ref<VectorXd> omega,
                                const Ref<const MatrixXd>& G,
                                const Ref<const VectorXd>& delta,
                                const Ref<const VectorXd>& epsilon) {
    // Need to have a valid starting value if the last one is not valid in MCMC
    if (omega != omega) { // if there is nan in omega
      // std::cout << "EvalOmegas: resetting omega_." << std::endl;
      omega = omega_init_;
    }
    int em_iter;
    double em_err;
    VectorXd weights(n_obs2_);
    eval_weights(weights, delta, epsilon, omega);
    double sum_weights = weights.head(n_obs2_).sum();
    norm_weights_.head(n_obs2_) /= sum_weights;
    double logel_old = GEL.logel_omega(omega, norm_weights_.head(n_obs2_), sum_weights);
    double logel = logel_old;
    int ii;
    for(ii=0; ii<max_iter_em_; ii++) {
      // E-step:
      eval_weights(weights, delta, epsilon, omega); // assigns weights according to epsilons
      sum_weights = weights.head(n_obs2_).sum();
      norm_weights_.head(n_obs2_) /= sum_weights;
      // M-step:
      VectorXd lambda;
      GEL.lambda_nr(lambda, G, norm_weights_.head(n_obs2_));
      // GEL.get_diag(nr_iter, nr_err);
      GEL.omega_hat(omega, lambda, G, norm_weights_.head(n_obs2_));
      logel = GEL.logel_omega(omega, norm_weights_.head(n_obs2_), sum_weights);
      em_err = abs(logel-logel_old); // absolute error in log EL
      if (em_err < abs_tol_) break;
      logel_old = logel;
    }
    em_iter = ii;
    if (em_iter == max_iter_em_ && em_err > abs_tol_) {
      // TODO: maybe should assign nan elsewhere
      for (int ii=0; ii<n_obs2_; ii++) {
        omega(ii) = std::numeric_limits<double>::quiet_NaN();
      }
    }
    return;
  }
  
  inline double CensEL::logel_omega(const Ref<const VectorXd>& delta,
                                    const Ref<const VectorXd>& epsilon,
                                    const Ref<const VectorXd>& omega) {
    VectorXd weights(n_obs2_);
    eval_weights(weights, delta, epsilon, omega);
    double sum_weights = weights.head(n_obs2_).sum();
    norm_weights_.head(n_obs2_) /= sum_weights;
    double logel = GEL.logel_omega(omega, norm_weights_, sum_weights);
    return logel;
  }

}

#endif