/// @file gen_el.h

#ifndef GEN_EL_H
#define GEN_EL_H

#include <RcppEigen.h>
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
    // for support modification
    int n_obs1_; // n_obs1_ = n_obs_+1: for memory allocation
    int n_obs2_; // n_obs2_ = n_obs_+support: used n_obs_ in calculation
    // logstar parameters
    VectorXd tstar_;
    VectorXd astar_;
    VectorXd bstar_;
    // VectorXd cstar_;
    VectorXd norm_weights_;
    // double trunc_, aa_, bb_, cc_; 
    // internal variables for Newton-Raphson 
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
    bool supp_adj_; // whether do support correction or not
    double supp_a_; // tuning parameter for support correction
    double weight_adj_; // weight of the additional observation under support correction
    VectorXd lambda0_; // initial lambda in Newton-Raphson iterations

    // functions to perform support correction.
    Ref<const MatrixXd> supp_G(const Ref<const MatrixXd>& G);
    Ref<const VectorXd> supp_norm_weights(const Ref<const VectorXd>& norm_weights);
    // implementation versions of functions to avoid defining pointer to matrix
    void omega_hat_impl(Ref<VectorXd> omega,
                        const Ref<const VectorXd>& lambda,
                        const Ref<const MatrixXd>& G,
                        const Ref<const VectorXd>& norm_weights);
    void lambda_nr_impl(Ref<VectorXd> lambda,
                        const Ref<const MatrixXd>& G,
                        const Ref<const VectorXd>& norm_weights);
    double logel_impl(const Ref<const MatrixXd>& G,
                      const Ref<const VectorXd>& norm_weights,
                      double sum_weights);
      
  public:

    /// Constructor.
    GenEL(int n_obs, int n_eqs);

    /// Set the maximum number of Newton-Raphson iterations.
    void set_max_iter(int max_iter);
    /// Set the relative tolerance for the Newton-Raphson algorithm.
    void set_rel_tol(double rel_tol);
    /// Set the support adjustment flag and support adjustment parameters.
    void set_supp_adj(bool supp_adj, double a, double weight_adj);
    /// Set the support adjustment flag.
    void set_supp_adj(bool supp_adj);
    /// Set the support adjustment tuning parameter.
    void set_supp_adj_a(double a);
    /// Set the weight for the additional fake observation under support adjustment.
    void set_weight_adj(double weight_adj);
    /// Set the initial value for Newton-Raphson algorithm.
    void set_lambda0(const Ref<const VectorXd>& lambda0);
    /// Get the diagnostics for the last Newton-Raphson run.
    void get_diag(int& nr_iter, double& nr_err);
    /// Get support adjustment flag.
    bool get_supp_adj();
    /// Get number of observations.
    int get_n_obs();
    /// Get number of estimating equations.
    int get_n_eqs();
    /// Get the maximum number of iterations.
    int get_max_iter();
    /// Get the relative tolerance.
    double get_rel_tol();
    
    /// Solve the dual problem via Newton-Raphson algorithm.
    void lambda_nr(Ref<VectorXd> lambda, 
                   const Ref<const MatrixXd>& G,
                   const Ref<const VectorXd>& weights);
                   // const Ref<const VectorXd>& norm_weights);
    /// Check whether Newton-Raphson algorithm has converged.
    bool has_converged_nr();
    /// Calculate the profile probability weights.
    void omega_hat(Ref<VectorXd> omega,
                   const Ref<const VectorXd>& lambda,
                   const Ref<const MatrixXd>& G,
                   const Ref<const VectorXd>& weights);
                   // const Ref<const VectorXd>& norm_weights);
    /// Calculate the empirical loglikelihood given the probability weights.
    double logel_omega(const Ref<const VectorXd>& omega,
                       const Ref<const VectorXd>& weights,
                       // const Ref<const VectorXd>& norm_weights,
                       double sum_weights);
    /// Calculate the unweighted empirical loglikelihood.
    double logel(const Ref<const MatrixXd>& G);
    /// Calculate the weighted empirical loglikelihood.
    double logel(const Ref<const MatrixXd>& G,
                 const Ref<const VectorXd>& weights);
    /// Calculate the gradient of the unweighted empirical loglikelihood.
    double logel_grad(Ref<MatrixXd> dldG, 
                      const Ref<const MatrixXd>& G);
    /// Calculate the gradient of the weighted empirical loglikelihood.
    double logel_grad(Ref<MatrixXd> dldG, 
                      const Ref<const MatrixXd>& G,
                      const Ref<const VectorXd>& weights);
    /// Calculate the gradient of the empirical loglikelihood.
    void logel_grad(Ref<MatrixXd> dldG,
                    const Ref<const VectorXd>& omega,
                    const Ref<const VectorXd>& lambda,
                    double sum_weights);
  };

  /// @param[in] n_obs    Number of observations.
  /// @param[in] n_eqs    Number of estimating equations.
  inline GenEL::GenEL(int n_obs, int n_eqs) {
    // assign internal values
    n_obs_ = n_obs;
    n_eqs_ = n_eqs;
    n_obs1_ = n_obs_+1; // space for augmented G
    tstar_ = VectorXd::Zero(n_obs1_);
    astar_ = VectorXd::Zero(n_obs1_);
    bstar_ = VectorXd::Zero(n_obs1_);
    norm_weights_ = VectorXd::Zero(n_obs1_);
    // memory allocation
    omega_ = VectorXd::Constant(n_obs1_, 1.0/(double)n_obs1_); // Initialize to 1/n_obs_
    Gaug_ = MatrixXd::Zero(n_eqs_,n_obs1_); // augmented G
    GGt_= MatrixXd::Zero(n_eqs_, n_eqs_); // columnwise outer product
    lambda_ = VectorXd::Zero(n_eqs_); // Initialize to all 0's
    lambda_new_ = VectorXd::Zero(n_eqs_);
    Q1_ = VectorXd::Zero(n_eqs_);
    Q2_ = MatrixXd::Zero(n_eqs_,n_eqs_);
    Glambda_ = VectorXd::Zero(n_obs1_);
    rho_ = VectorXd::Zero(n_obs1_);
    Q2llt_.compute(MatrixXd::Identity(n_eqs_,n_eqs_));
    // initialize NR options
    nr_iter_ = 0;
    nr_err_ = 0.0;
    set_max_iter(100);
    set_rel_tol(1e-7);
    set_supp_adj(false);
    set_weight_adj(1.0);
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
  /// @param[in] supp_adj   Whether or not to enable support adjustment.
  /// @param[in] a          Support adjustment factor.  Defaults to `max(1.0, log(n_obs)/2)`.
  /// @param[in] weight_adj Weight for the additional observation when calculating weighted log EL.
  inline void GenEL::set_supp_adj(bool supp_adj, double a, double weight_adj) {
    supp_adj_ = supp_adj;
    supp_a_ = a;
    weight_adj_ = weight_adj;
    n_obs2_ = n_obs_+supp_adj_;
    return;
  }
  
  /// Reference: J. Chen, A. M. Variyath, and B. Abraham. Adjusted empirical likelihood and its properties. Journal of Computational and Graphical Statistics, 17(2):426–443, 2008).
  ///
  /// @param[in] supp_adj   Whether or not to enable support adjustment.
  inline void GenEL::set_supp_adj(bool supp_adj) {
    set_supp_adj(supp_adj, std::max(1.0,0.5*log(n_obs_)), 1.0);
    return;
  }
  
  /// Reference: J. Chen, A. M. Variyath, and B. Abraham. Adjusted empirical likelihood and its properties. Journal of Computational and Graphical Statistics, 17(2):426–443, 2008).
  ///
  /// @param[in] a          Support adjustment factor.  Defaults to `max(1.0, log(n_obs)/2)`.
  inline void GenEL::set_supp_adj_a(double a) {
    supp_a_ = a;
    return;
  }
  
  /// @param[in] weight_adj Weight for the additional observation when calculating weighted log EL.
  inline void GenEL::set_weight_adj(double weight_adj) {
    weight_adj_ = weight_adj;
    return;
  }

  /// @param[in] lambda0 Initial value of lambda of length `n_eqs`.
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

  /// @return Whether the desired error tolerance `rel_tol` has been reached in less than `max_iter` Newton-Raphson steps.
  /// @warning Assumes that lambda_nr has been run at least once. 
  inline bool GenEL::has_converged_nr() {
    return !((nr_err_ > rel_tol_) && (nr_iter_ == max_iter_));
  }

  /// @param[in] G Moment matrix of size `n_eqs x n_obs` or `n_eqs x (n_obs + 1)`. If `supp_adj = false`, the former is required.  If `supp_adj = true` and the former is provided, support adjustment is performed.  If `supp_adj = true` and `G.cols() == n_obs + 1`, assumes that support has already been corrected.
  inline Ref<const MatrixXd> GenEL::supp_G(const Ref<const MatrixXd>& G) {
    if(supp_adj_ && G.cols() == n_obs_) {
      Gaug_.leftCols(n_obs_) = G;
      adj_G(Gaug_, supp_a_);
      return Ref<const MatrixXd>(Gaug_);
    } else {
      return Ref<const MatrixXd>(G);
    }
  }

  // /// @param[in] norm_weights Vector of normalized weights (sum to 1) of length `n_obs` or `n_obs + 1`. 
  // ///  If `supp_adj = false`, the former is required.  
  // ///  If `supp_adj = true` and the former is provided, weight renormalization is performed. i.e., 
  // ///  The the first `n_obs` entries are renormalized to sum to `n_obs/(n_obs+1.0)` and 
  // ///  the last entry will have a value `1.0/(n_obs+1.0)`.
  // ///  If `supp_adj = true` and `norm_weights.size() == n_obs + 1`, assumes that 
  // ///  weights have already been renormalized. 
  // ///  The return value still sums to one.
  // inline Ref<const VectorXd> GenEL::supp_norm_weights(const Ref<const VectorXd>& norm_weights) {
  //   if(supp_adj_ && norm_weights.size() == n_obs_) {
  //     norm_weights_.head(n_obs_) = n_obs_/(n_obs_+1.0) * norm_weights;
  //     norm_weights_(n_obs_) = 1.0/(n_obs_+1.0);
  //     return Ref<const VectorXd>(norm_weights_);
  //   } else {
  //     return Ref<const VectorXd>(norm_weights);
  //   }
  // }

  inline Ref<const VectorXd> GenEL::supp_norm_weights(const Ref<const VectorXd>& weights) {
    if(supp_adj_ && weights.size() == n_obs_) {
      norm_weights_.head(n_obs_) = weights;
      norm_weights_(n_obs_) = weight_adj_;
      // std::cout << "norm_weights_ before normalization = " << norm_weights_.transpose() << std::endl;
      norm_weights_ = norm_weights_/norm_weights_.sum();
      return Ref<const VectorXd>(norm_weights_);
    } else {
      // std::cout << "supp_norm_weights: only normalize." << std::endl;
      // std::cout << "weights = " << weights.transpose() << std::endl;
      norm_weights_.head(weights.size()) = weights/weights.sum();
      return Ref<const VectorXd>(norm_weights_.head(weights.size()));
    }
  }
  
  inline void GenEL::lambda_nr_impl(Ref<VectorXd> lambda,
                                    const Ref<const MatrixXd>& G,
                                    const Ref<const VectorXd>& norm_weights) {
    int n_obs = G.cols();
    lambda = lambda0_; // set to initial value
    // logstar constants
    // assume the weights are normalized to sum to one
    bstar_ = norm_weights.array().inverse();
    astar_ = -.5 * bstar_.array().square();
    bstar_ *= 2.0;
    // newton-raphson loop
    int ii, jj;
    for(ii=0; ii<max_iter_; ii++) {
      // Q1 and Q2
      Glambda_.noalias() = lambda.transpose() * G;
      Glambda_= 1.0 - Glambda_.array();
      Q2_.fill(0.0);
      for(jj=0; jj<n_obs; jj++) {
        rho_(jj) = log_star1(Glambda_(jj),
             norm_weights(jj), astar_(jj), bstar_(jj));
        GGt_.noalias() = G.col(jj) * G.col(jj).transpose();
        Q2_ -= norm_weights(jj) *
          log_star2(Glambda_(jj), norm_weights(jj), astar_(jj)) * GGt_;
      }
      // update lambda
      rho_.array() *= norm_weights.array();
      Q1_.noalias() = G * rho_;
      Q2llt_.compute(Q2_);
      Q2llt_.solveInPlace(Q1_);
      lambda_new_ = lambda - Q1_;
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
  /// @param[in] G Moment matrix of size `n_eqs x n_obs` or `n_eqs x (n_obs + supp_adj)`.  If `supp_adj = false`, the former is required.  If `supp_adj = true` and the former is provided, support adjustment is performed.  If `supp_adj = true` and `G.cols() == n_obs + 1`, assumes that support has already been corrected. 
  /// @param[in] norm_weights Vector of weights of length `n_obs` or `n_obs + 1` following same logic as above.  Normalized to sum to one.
  ///
  /// @note Due to the difficulty of creating Eigen "pointers", current hack is to explicitly tell lambda_nr_impl whether to use external or internal `G` and `norm_weights`.
  inline void GenEL::lambda_nr(Ref<VectorXd> lambda,
                               const Ref<const MatrixXd>& G,
                               const Ref<const VectorXd>& weights) {  
                               // const Ref<const VectorXd>& norm_weights) {  
    Ref<const MatrixXd> G_eff = supp_G(G);
    Ref<const VectorXd> norm_weights_eff = supp_norm_weights(weights);
    // std::cout << "norm_weights_eff = " << norm_weights_eff.transpose() << std::endl;
    lambda_nr_impl(lambda, G_eff, norm_weights_eff);
    return;
  }

  inline void GenEL::omega_hat_impl(Ref<VectorXd> omega,
                                    const Ref<const VectorXd>& lambda,
                                    const Ref<const MatrixXd>& G,
                                    const Ref<const VectorXd>& norm_weights) {
    Glambda_.noalias() = lambda.transpose() * G;
    omega = (1.0-Glambda_.array()).inverse() * norm_weights.array();
    // note: when lambda is the optimal value denominator is n_obs.
    // however things get more mixed up with support correction,
    // so explicitly normalize for safety.
    omega.array() /= omega.array().sum();
    return;
  }

  
  /// @param[out] omega Probability vector of length `n_obs + supp_adj`.
  /// @param[in] lambda Dual problem vector of size `n_eqs`.  
  /// @param[in] G Moment matrix of size `n_eqs x n_obs` or `n_eqs x (n_obs + supp_adj)`.  If `supp_adj = false`, the former is required.  If `supp_adj = true` and the former is provided, support adjustment is performed.  If `supp_adj = true` and `G.cols() == n_obs + 1`, assumes that support has already been corrected. 
  /// @param[in] norm_weights Vector of weights of length `n_obs` or `n_obs + supp_adj` following same logic as above.  Normalized to sum to one.
  inline void GenEL::omega_hat(Ref<VectorXd> omega,
                               const Ref<const VectorXd>& lambda,
                               const Ref<const MatrixXd>& G,
                               const Ref<const VectorXd>& weights) {
                               // const Ref<const VectorXd>& norm_weights) {
    Ref<const MatrixXd> G_eff = supp_G(G);
    Ref<const VectorXd> norm_weights_eff = supp_norm_weights(weights);
    omega_hat_impl(omega, lambda, G_eff, norm_weights_eff); 
    return;
  }

  /// @param[in] omega Probability vector of length `n_obs + supp_adj`.
  /// @param[in] norm_weights Vector of weights of length `n_obs` or `n_obs + supp_adj` following same logic as above.  Normalized to sum to one.
  /// @param[in] sum_weights Sum of the weights prior to normalization.
  /// @return Value of empirical loglikelihood, which is `sum(log(omega))`.
  inline double GenEL::logel_omega(const Ref<const VectorXd>& omega,
                                   // const Ref<const VectorXd>& norm_weights,
                                   const Ref<const VectorXd>& weights,
                                   double sum_weights) {
    Ref<const VectorXd> norm_weights_eff = supp_norm_weights(weights);
    // std::cout << "sum_weights = " << sum_weights << std::endl;
    // std::cout << "omega = " << omega.transpose() << std::endl;
    // std::cout << "norm_weights_eff = " << norm_weights_eff.transpose() << std::endl;
    return sum_weights * (omega.array().log() * norm_weights_eff.array()).sum();
  }

  inline double GenEL::logel_impl(const Ref<const MatrixXd>& G,
                                  const Ref<const VectorXd>& norm_weights,
                                  double sum_weights) {
    lambda_nr(lambda_, G, norm_weights);
    omega_hat(omega_.head(n_obs2_), lambda_, G, norm_weights);
    return logel_omega(omega_.head(n_obs2_), norm_weights, sum_weights);
  }

  /// @param[in] G Moment matrix of size `n_eqs x n_obs`.
  /// @param[in] weights Vector of weights of length `n_obs`.  Not assumed to be normalized.
  /// @return Value of empirical loglikelihood.
  inline double GenEL::logel(const Ref<const MatrixXd>& G,
                             const Ref<const VectorXd>& weights) {
    Ref<const MatrixXd> G_eff = supp_G(G);
    double sum_weights = weights.sum();
    norm_weights_.head(n_obs_) = weights;
    if (supp_adj_) {
      norm_weights_(n_obs_) = weight_adj_;
      sum_weights += weight_adj_;
    }
    norm_weights_.head(n_obs2_) /= sum_weights;
    return logel_impl(G_eff, norm_weights_.head(n_obs2_), sum_weights);
  }

  /// @param[in] G Moment matrix of size `n_eqs x n_obs`.
  /// @return The value of the empirical loglikelihood.
  inline double GenEL::logel(const Ref<const MatrixXd>& G) {
    Ref<const MatrixXd> G_eff = supp_G(G);
    double sum_weights = double(n_obs2_);
    norm_weights_.head(n_obs2_) = VectorXd::Ones(n_obs2_)/sum_weights;
    return logel_impl(G_eff, norm_weights_.head(n_obs2_), sum_weights);
  }

  /// @param[out] dldG Matrix of size `n_eqs x n_obs` containing the gradient of logel with respect to G.
  /// @param[in] omega Probability vector of length `n_obs + supp_adj`.
  /// @param[in] lambda Dual problem vector of size `n_eqs`.
  /// @param[in] sum_weights Sum of the weights prior to normalization.
  inline void GenEL::logel_grad(Ref<MatrixXd> dldG,
                                const Ref<const VectorXd>& omega,
                                const Ref<const VectorXd>& lambda,
                                double sum_weights) {
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
  inline double GenEL::logel_grad(Ref<MatrixXd> dldG,
                                  const Ref<const MatrixXd>& G,
                                  const Ref<const VectorXd>& weights) {
    double sum_weights = weights.sum();
    if (supp_adj_) {
      sum_weights += weight_adj_;
    }
    double ll = logel(G, weights);
    // use internal values of omega_ and lambda_
    logel_grad(dldG, omega_.head(n_obs2_), lambda_, sum_weights);
    return ll;
  }

  /// @param[out] dldG Matrix of size `n_eqs x n_obs` containing the gradient of logel with respect to G.
  /// @param[in] G Moment matrix of size `n_eqs x n_obs`.
  /// @return For convenience, returns the value of the empirical loglikelihood.
  inline double GenEL::logel_grad(Ref<MatrixXd> dldG,
                                  const Ref<const MatrixXd>& G) {
    double ll = logel(G);
    // use internal values of omega_ and lambda_
    logel_grad(dldG, omega_.head(n_obs2_), lambda_, double(n_obs2_));
    return ll;
  }
  
  inline bool GenEL::get_supp_adj() {
    return supp_adj_;
  }
  
  inline int GenEL::get_n_obs() {
    return n_obs_;
  }
  
  inline int GenEL::get_n_eqs() {
    return n_eqs_;
  }
  
  inline int GenEL::get_max_iter() {
    return max_iter_;
  }
  
  inline double GenEL::get_rel_tol() {
    return rel_tol_;
  }
  
} // end namespace flexEL

#endif
