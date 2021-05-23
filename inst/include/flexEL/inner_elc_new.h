/// @file inner_elc.h

#ifndef INNER_ELC_H
#define INNER_ELC_H

#include <RcppEigen.h>
#include "inner_el.h"
#include "sort_order.h"
#include "ind_smooth.h" 

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
    // for continuity correction
    bool smooth_;
    double smooth_s_;
    
    // internal variables for EM algorithm
    VectorXd omega_init_; // initial value for omegas, in case of reset
    VectorXd delta_; // censoring indicators
    VectorXd epsilon_;
    VectorXi eps_ord_; 
    VectorXd norm_weights_;
    VectorXd psot_; // partial sum of omegatildas
    VectorXd psos_vec_; // partial sum of omegas up to each omega in the order of epsilon
    int em_iter_;
    double em_err_;

    // internals for EM options
    int max_iter_em_;
    double abs_tol_; // absolute tolerance
    
    // helper function for eval_weights: calculate partial sum of omegas
    void eval_pso(double &psos,
                  const int ii,
                  const Ref<const VectorXd>& omega); // Partial sum of omega_jj s.t. eps_jj >= eps_ii
    // helper function for eval_weights: calculate partial sum of omegas
    void eval_pso_smooth(double &psos,
                         const int ii,
                         const Ref<const VectorXd>& omega,
                         const Ref<const VectorXd>& epsilon);
    
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
    /// Set the discontinuity correction flag
    void set_smooth(bool smooth, double s);
    void set_smooth(bool smooth);
    /// Set the initial value for Newton-Raphson algorithm.
    void set_lambda0(const Ref<const VectorXd>& lambda0);
    /// Get the diagnostics for the EM run.
    void get_diag(int& em_iter, double& em_err);
    /// Get support adjustment flag.
    bool get_supp_adj();
    // Get continuity correction flag.
    bool get_smooth();
    /// Get number of observations.
    int get_n_obs();
    /// Get number of estimating equations.
    int get_n_eqs();
    /// Get the maximum number of iterations for Newton-Raphson.
    int get_max_iter_nr();
    /// Get the relative tolerance.
    double get_rel_tol();
    /// Get the maximum number of iterations for EM.
    int get_max_iter_em();
    /// Get the relative tolerance.
    double get_abs_tol();
    
    /// Calculate weights for EM algorithm.
    void eval_weights(Ref<VectorXd> weights,
                      const Ref<const VectorXd>& delta,
                      const Ref<const VectorXd>& epsilon,
                      const Ref<const VectorXd>& omega);
    
    /// Calculate weights for EM algorithm with continuity correction.
    void eval_weights_smooth(Ref<VectorXd> weights,
                             const Ref<const VectorXd>& delta,
                             const Ref<const VectorXd>& epsilon,
                             const Ref<const VectorXd>& omega);
    
    /// Calculate the profile probability weights.
    void omega_hat(Ref<VectorXd> omega,
                   const Ref<const MatrixXd>& G,
                   const Ref<const VectorXd>& delta,
                   const Ref<const VectorXd>& epsilon);
    
    /// Calculate the empirical loglikelihood.
    double logel(const Ref<const MatrixXd>& G,
                 const Ref<const VectorXd>& delta,
                 const Ref<const VectorXd>& epsilon);
    
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
    smooth_ = false;
    smooth_s_ = 10;

    omega_init_ = VectorXd::Zero(n_obs1_).array() + 1.0/(double)n_obs1_; 
    delta_ = VectorXd::Zero(n_obs1_);
    epsilon_ = VectorXd::Zero(n_obs1_);
    eps_ord_ = VectorXi::Zero(n_obs1_);
    norm_weights_ = VectorXd::Zero(n_obs1_);
    psot_ = VectorXd::Zero(n_obs1_);
    psos_vec_ = VectorXd::Zero(n_obs1_);
    
    set_max_iter_em(100);
    set_abs_tol(1e-3);
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
  
  /// @param[in] abs_tol Relative tolerance for the EM algorithm.
  inline void CensEL::set_abs_tol(double abs_tol) {
    abs_tol_ = abs_tol;
    return;
  }
  
  /// @param[in] supp_adj Whether do support correction or not.
  /// @param[in] a        A multiplier to the additional estimating equation.
  inline void CensEL::set_supp_adj(bool supp_adj, double a) {
    supp_adj_ = supp_adj; 
    supp_a_ = a; 
    n_obs2_ = n_obs_+supp_adj_;
    GEL.set_supp_adj(supp_adj_, supp_a_);
    if (supp_adj) {
      epsilon_.tail(1)(0) = -INFINITY;
      delta_.tail(1)(0) = 0;
    }
    return;
  }
  
  /// @param[in] supp_adj Whether do support correction or not.
  inline void CensEL::set_supp_adj(bool supp_adj) {
    set_supp_adj(supp_adj, std::max(1.0,0.5*log(n_obs_)));
    return;
  }
  
  /// @param[in] smooth Whether do continuity correction or not.
  /// @param[in] s      Tuning parameter for continuity correction.
  inline void CensEL::set_smooth(bool smooth, double s) {
    smooth_ = smooth; 
    smooth_s_ = s; 
    return;
  }
  
  /// @param[in] smooth Whether do continuity correction or not.
  inline void CensEL::set_smooth(bool smooth) {
    set_smooth(smooth, 10); 
    return;
  }
  
  /// @param[in] lambda0 Initial value of lambda of length `n_eqs`.
  inline void CensEL::set_lambda0(const Ref<const VectorXd>& lambda0) {
    GEL.set_lambda0(lambda0);
    return;
  }
  
  /// @param[out] em_iter Number of EM iterations.
  /// @param[out] em_err Maximum relative difference between the value of logel in the last two EM steps.
  inline void CensEL::get_diag(int& em_iter, double& em_err) {
    // GEL.get_diag(nr_iter, nr_err);
    em_iter = em_iter_;
    em_err = em_err_;
    return;
  }
  
  /// @param[out] psos Partial sum of omega.
  /// @param[in] ii Index for ii-th largest epsilon.
  /// @param[in] omega Probability vector of length `n_obs + supp_adj`.
  inline void CensEL::eval_pso(double &psos,
                               const int ii,
                               const Ref<const VectorXd>& omega) {
    psos = 0;
    int kk;
    for (int jj=n_obs2_-1; jj>=0; jj--) {
      // Note: eps_ord_ corresponds to epsilons in ascending order
      kk = eps_ord_(jj); // kk is the index of the jj-th largest epsilon
      psos += omega(kk);
      if (kk == ii) break; // until (backwardly) reaching ii-th largest epsilon 
    }
    return;
  }
  
  inline void CensEL::eval_pso_smooth(double &psos,
                                      const int ii,
                                      const Ref<const VectorXd>& omega,
                                      const Ref<const VectorXd>& epsilon) {
    // std::cout << "epsilon = " << epsilon.transpose() << std::endl;
    if (supp_adj_ && ii == (n_obs2_-1)) {
      psos = omega.head(n_obs2_-1).sum() + 0.5*omega(n_obs2_-1);
    }
    else {
      for (int jj=0; jj<n_obs2_; jj++) {
        // std::cout << "smooth_s_ = " << smooth_s_ << std::endl;
        psos += ind_smooth(epsilon(ii)-epsilon(jj), smooth_s_)*omega(jj);
      }
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
    double psos; // must init to 0
    delta_.head(n_obs_) = delta; 
    epsilon_.head(n_obs_) = epsilon;
    // std::cout << "n_obs2_ = " << n_obs2_ << std::endl;
    eps_ord_.head(n_obs2_) = sort_inds(epsilon_.head(n_obs2_)); 
    // std::cout << "eps_ord_.head(n_obs2_) = " << eps_ord_.head(n_obs2_).transpose() << std::endl;
    // std::cout << "eps_ord_ = " << eps_ord_.transpose() << std::endl;
    for (int ii=0; ii<n_obs2_; ii++) {
      for (int jj=0; jj<n_obs2_; jj++) {
        kk = eps_ord_(jj);
        if (delta_(kk) == 0) {
          eval_pso(psos, kk, omega);
          // std::cout << "psos = " << psos << std::endl;
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
  
  inline void CensEL::eval_weights_smooth(Ref<VectorXd> weights,
                                          const Ref<const VectorXd>& delta,
                                          const Ref<const VectorXd>& epsilon,
                                          const Ref<const VectorXd>& omega) {
    psot_.fill(0.0); 
    psos_vec_.fill(0.0);
    delta_.head(n_obs_) = delta; 
    epsilon_.head(n_obs_) = epsilon; 
    // std::cout << "epsilon_.head(n_obs2_) = " << epsilon_.head(n_obs2_).transpose() << std::endl;
    for (int ii=0; ii<n_obs2_; ii++) {
      eval_pso_smooth(psos_vec_(ii), ii, omega, epsilon_.head(n_obs2_));
    }
    // std::cout << "psos_vec_ = " << psos_vec_.transpose() << std::endl;

    if (supp_adj_) {
      for (int jj=0; jj<n_obs2_; jj++) {
        // std::cout << "jj = " << jj << std::endl;
        for (int kk=0; kk<n_obs2_; kk++) {
          // std::cout << "  kk = " << kk << std::endl;
          if (jj == n_obs2_-1 && kk == n_obs2_-1) {
            psot_(jj) += (1-delta_(kk))*ind_smooth(0.0,smooth_s_)*omega(jj)/psos_vec_(kk);
          }
          else {
            psot_(jj) += (1-delta_(kk))*ind_smooth(epsilon_(kk)-epsilon_(jj),smooth_s_)*omega(jj)/psos_vec_(kk);
          }
        }
      }
    }
    else {
      for (int jj=0; jj<n_obs2_; jj++) {
        for (int kk=0; kk<n_obs2_; kk++) {
          psot_(jj) += (1-delta_(kk))*ind_smooth(epsilon_(kk)-epsilon_(jj),smooth_s_)*omega(jj)/psos_vec_(kk);
        }
      }
    }
    // std::cout << "psot_ = " << psot_.transpose() << std::endl;
    weights = delta_.array() + psot_.array();
  }
  
  /// @param[out] omega  Probability vector of length `n_obs + supp_adj`.
  /// @param[in] G       Moment matrix of size `n_eqs x n_obs` or `n_eqs x (n_obs + supp_adj)`.  If `supp_adj = false`, the former is required.  If `supp_adj = true` and the former is provided, support adjustment is performed.  If `supp_adj = true` and `G.cols() == n_obs + 1`, assumes that support has already been corrected. 
  /// @param[in] delta   Vector of censoring indicators of length `n_obs` or `n_obs + supp_adj` following same logic as above.
  inline void CensEL::omega_hat(Ref<VectorXd> omega,
                                const Ref<const MatrixXd>& G,
                                const Ref<const VectorXd>& delta,
                                const Ref<const VectorXd>& epsilon) {
    // std::cout << "omega_hat: omega = " << omega.transpose() << std::endl;
    // Need to have a valid starting value if the last one is not valid in MCMC
    if (omega != omega) { // if there is nan in omega
      // std::cout << "omega_hat: resetting omega_." << std::endl;
      omega = omega_init_;
    }
    int em_iter;
    double em_err;
    VectorXd weights(n_obs2_);
    if (!smooth_) {
      eval_weights(weights, delta, epsilon, omega);
      // std::cout << "delta = " << omega.transpose() << std::endl;
      // std::cout << "epsilon = " << omega.transpose() << std::endl;
      // std::cout << "initial omega = " << omega.transpose() << std::endl;
      // std::cout << "initial weights = " << weights.transpose() << std::endl;
    } else{
      eval_weights_smooth(weights, delta, epsilon, omega);
    }
    double sum_weights = weights.head(n_obs2_).sum();
    norm_weights_.head(n_obs2_) = weights/sum_weights;
    // std::cout << "initial norm_weights_.head(n_obs2_) = " << norm_weights_.head(n_obs2_).transpose() << std::endl;
    double logel_old = GEL.logel_omega(omega, norm_weights_.head(n_obs2_), sum_weights);
    double logel = logel_old;
    // std::cout << "initial logel = " << logel << std::endl;
    int ii;
    VectorXd lambda(n_eqs_);
    // std::cout << "before EM loop" << std::endl;
    // std::cout << "norm_weights_ = " << norm_weights_.transpose() << std::endl;
    for(ii=0; ii<max_iter_em_; ii++) {
      // std::cout << "EM ii = " << ii << std::endl;
      // E-step:
      if (!smooth_) {
        eval_weights(weights, delta, epsilon, omega);
      } else{
        eval_weights_smooth(weights, delta, epsilon, omega);
      }
      sum_weights = weights.head(n_obs2_).sum();
      norm_weights_.head(n_obs2_) = weights/sum_weights;
      // std::cout << "norm_weights_ = " << norm_weights_.transpose() << std::endl;
      // M-step:
      GEL.lambda_nr(lambda, G, norm_weights_.head(n_obs2_));
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
    em_iter_ = em_iter;
    em_err_ = em_err;
    return;
  }
  
  // TODO: this doesn't give correct result..
  // /// @param[in] delta Censoring indicator vector of length `n_obs + supp_adj`.
  // /// @param[in] epsilon Residual vector of length `n_obs + supp_adj`.
  // /// @param[in] omega Probability vector of length `n_obs + supp_adj`.
  // inline double CensEL::logel_omega(const Ref<const VectorXd>& delta,
  //                                   const Ref<const VectorXd>& epsilon,
  //                                   const Ref<const VectorXd>& omega) {
  //   VectorXd weights(n_obs2_);
  //   // eval_weights(weights, delta, epsilon, omega);
  //   if (!smooth_) {
  //     eval_weights(weights, delta.head(n_obs_), epsilon.head(n_obs_), omega);
  //     std::cout << "weights = " << weights.transpose() << std::endl;
  //   } else{
  //     eval_weights_smooth(weights, delta.head(n_obs_), epsilon.head(n_obs_), omega);
  //   }
  //   double sum_weights = weights.head(n_obs2_).sum();
  //   norm_weights_.head(n_obs2_) /= sum_weights;
  //   double logel = GEL.logel_omega(omega, norm_weights_, sum_weights);
  //   return logel;
  // }
  
  /// @param[in] G       Moment matrix of size `n_eqs x n_obs`. 
  /// @param[in] delta   Censoring indicator vector of length `n_obs`.
  /// @param[in] epsilon Residual vector of length `n_obs`.
  inline double CensEL::logel(const Ref<const MatrixXd>& G,
                              const Ref<const VectorXd>& delta,
                              const Ref<const VectorXd>& epsilon) {
    // TODO: dimension of G works?
    VectorXd omega = Eigen::VectorXd::Constant(n_obs2_, 1.0/(n_obs2_));
    delta_.head(n_obs_) = delta;
    epsilon_.head(n_obs_) = epsilon;
    omega_hat(omega, G, delta_.head(n_obs2_), epsilon_.head(n_obs2_));
    
    VectorXd weights(n_obs2_);
    if (!smooth_) {
      eval_weights(weights, delta.head(n_obs_), epsilon.head(n_obs_), omega);
      // std::cout << "weights = " << weights.transpose() << std::endl;
    } else{
      eval_weights_smooth(weights, delta.head(n_obs_), epsilon.head(n_obs_), omega);
    }
    if (supp_adj_) GEL.set_weight_adj(weights(n_obs_));
    double logel = GEL.logel(G, weights.head(n_obs_));
    // omega_hat(omega, G, delta_.head(n_obs2_), epsilon_.head(n_obs2_));
    // std::cout << "omega = " << omega.transpose() << std::endl;
    // std::cout << "delta_.head(n_obs2_) = " << delta_.head(n_obs2_).transpose() << std::endl;
    // std::cout << "epsilon_.head(n_obs2_) = " << epsilon_.head(n_obs2_).transpose() << std::endl;
    // double logel = logel_omega(delta_.head(n_obs2_), epsilon_.head(n_obs2_), omega);
    return logel;
  }
  
  inline bool CensEL::get_supp_adj() {
    return supp_adj_;
  }
  
  inline bool CensEL::get_smooth() {
    return smooth_;
  }
  
  inline int CensEL::get_n_obs() {
    return n_obs_;
  }
  
  inline int CensEL::get_n_eqs() {
    return n_eqs_;
  }
  
  inline int CensEL::get_max_iter_nr() {
    return GEL.get_max_iter();
  }
  
  inline double CensEL::get_rel_tol() {
    return GEL.get_rel_tol();
  }
  
  inline int CensEL::get_max_iter_em() {
    return max_iter_em_;
  }
  
  inline double CensEL::get_abs_tol() {
    return abs_tol_;
  }

}

#endif
