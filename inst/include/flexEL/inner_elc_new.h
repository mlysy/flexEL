/// @file inner_elc.h

#ifndef INNER_ELC_H
#define INNER_ELC_H

#include <RcppEigen.h>
#include "adj_G.h" // for support adjustment
#include "log_star.h"

namespace flexEL {

  using namespace Eigen;

  class CensEL {
    
  private:
    
    GenEL GEL;
    
    // internal variables for EM algorithm
    VectorXd delta_; // censoring indicators
    VectorXd weights_; // weights in weighted maximum log EL
    VectorXd omega_; // empirical distribution
    VectorXd omega_init_; // initial value for omegas, in case of reset
    VectorXd omegaps_; // partial sum of omegas
    VectorXd epsilon_; // residuals used for ordering in EM
    VectorXi eps_ord_; // order of epsilons
    VectorXd psot_; // partial sum of omegatildas
    VectorXd pso_; // partial sum of omegas
    ArrayXd logel_arr_; // an array to assist logel calculation 
    
    // internals for EM options
    double abs_tol_; // absolute tolerance
    
    // helper function for eval_weights: calculate partial sum of omegas
    double eval_pso(const int ii); // Partial sum of omegas_jj s.t. eps_jj >= eps_ii
    
    
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
    
    /// Calculate weights according to epsilons
    void eval_weights();
    
  };
  
  /// @param[in] n_obs    Number of observations.
  /// @param[in] n_eqs    Number of estimating equations.
  inline CensEL::CensEL(int n_obs, int n_eqs) {
    
    GEL = GenEL(n_obs, n_eqs);
    
    // TODO: shall CensEL still have it's own n_eqs_, n_obs_, etc. or use them from GEL?
    
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
  
  /// @param[out] lambda Vector of length `n_eqs` containing the candidate solution to the dual optimization problem.
  /// @param[in] G Moment matrix of size `n_eqs x n_obs` or `n_eqs x (n_obs + supp_adj)`.  If `supp_adj = false`, the former is required.  If `supp_adj = true` and the former is provided, support adjustment is performed.  If `supp_adj = true` and `G.cols() == n_obs + 1`, assumes that support has already been corrected. 
  /// @param[in] norm_weights Vector of weights of length `n_obs` or `n_obs + 1` following same logic as above.  Normalized to sum to one.
  ///
  /// @note Due to the difficulty of creating Eigen "pointers", current hack is to explicitly tell lambda_nr_impl whether to use external or internal `G` and `norm_weights`.
  inline void CensEL::lambda_nr(Ref<VectorXd> lambda,
                               const Ref<const MatrixXd>& G,
                               const Ref<const VectorXd>& norm_weights) {  
    GEL.lambda_nr(lambda, G, norm_weights); // TODO: add methods to calculate weights for CensEL
    return;
  }
  
  // TODO: needs modification
  inline double CensEL::eval_pso(const int ii) {
    double psos = 0;
    int kk;
    for (int jj=n_obs2_-1; jj>=0; jj--) {
      // Note: eps_ord_ corresponds to epsilons in ascending order
      kk = eps_ord_(jj); // kk is the index of the jj-th largest epsilon
      psos += omega_(kk);
      if (kk == ii) break; // until (backwardly) reaching ii-th largest epsilon 
    }
    return psos;
  }
  
  // TODO: needs modification
  inline void CensEL::eval_weights() {
    // find the indices for increasing order of epsilons 
    psot_.fill(0.0);
    int kk;
    double psos;
    for (int ii=0; ii<n_obs2_; ii++) {
      for (int jj=0; jj<n_obs2_; jj++) {
        kk = eps_ord_(jj);
        if (delta_(kk) == 0) {
          psos = eval_pso(kk);
          // to prevent dividing by 0
          if (abs(psos) >= 1e-10) psot_(ii) += omegas_(ii)/psos;
          // else if (omegas_(ii) >= 1e-10 && EvalPSO(kk) < 1e-10) {
          //   // TODO: this means a problem
          //   // std::cout << "EvalWeights: dividing by 0 problem." << std::endl;
          // }
        }
        if (kk == ii) break;
      }
    }
    weights_.array() = delta_.array() + psot_.array();
  }

}

#endif