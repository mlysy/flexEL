/**
 * @file InnerEL.h
 * 
 * @brief Inner optimization for empirial likelihood problems.
 */

#ifndef INNEREL_h
#define INNEREL_h

#include <Rcpp.h>
#include <RcppEigen.h>
#include "BlockOuter.h" // columnwise outer product
#include "AdjG.h" // for support correction
// #include "MwgAdapt.h" // for adaptive mcmc

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

/* --------------------------------------------------------------------------- */

/**
 * @brief EL namespace
 * 
 * Wrap the exported library components into a namespace called \b el to avoid potential naming conflicts with other libraries or user-defined headers.
 */
namespace el {
  
  /**
   * @file       InnerEL.h
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
    // double MaxRelErr(const Ref<const VectorXd>& lambdaNew,
    //                  const Ref<const VectorXd>& lambdaOld);
    
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
    // void setTol(const int& max_iter, const double& rel_tol);
    
    // core caculations
    void LambdaNR(int& n_iter, double& max_iter); // Note: rel_tol and max_iter must be set before calling
    void EvalOmegas(); // empirical distribution
    double LogEL(); // log empirical likelihood
    
    // set and get functions 
    void set_lambda(const Ref<const VectorXd>& lambda); // assigned to lambdaNew
    void set_omegas(const Ref<const VectorXd>& omegas); 
    void set_G(const Ref<const MatrixXd>& G);
    VectorXd get_lambda(); 
    VectorXd get_omegas();
    MatrixXd get_G(); // TODO: is it better to return a reference?..
    Ref<MatrixXd> get_ref_G();
    
    // nBet, nGam and nQts are FOR THE MCMC SAMPELERS
    // using ELModel::nBet;
    // using ELModel::nGam;
    // using ELModel::nQts;
    // // posterior samplers:
    // void mwgStep(VectorXd &thetaCur, const int &idx, const double &mwgsd,
    //              bool &accept, double &LogELCur);
    // MatrixXd postSample(int nsamples, int nburn, MatrixXd ThetaInit,
    //                     const Ref<const MatrixXd>& MwgSds, 
    //                     MatrixXd &RvDoMcmc, MatrixXd &Paccept);
    // MatrixXd postSampleAdapt(int nsamples, int nburn, VectorXd thetaInit,
    //                          double *mwgSd, bool *rvDoMcmc, VectorXd &paccept);
  };

} // namespace el

/* --------------------------------------------------------------------------- */

// private functions

// maximum relative error in lambda
inline double el::InnerEL::MaxRelErr() {
  // TODO: added for numerical stability, what is a good tolerance to use ?
  if ((lambda_new_ - lambda_old_).array().abs().maxCoeff() < 1e-10) return(0);
  
  rel_err_ = ((lambda_new_ - lambda_old_).array() / (lambda_new_ + lambda_old_).array()).abs();
  return(rel_err_.maxCoeff());
}

// LogStar
inline double el::InnerEL::LogStar(double x) {
  if(x >= trunc_) {
    return(log(x));
  } else {
    return((aa_*x + bb_)*x + cc_);
  }
}

// d LogStar(x)/dx
inline double el::InnerEL::LogStar1(double x) {
  if(x >= trunc_) {
    return(1.0/x);
  } 
  else {
    // return(aa*x + bb);
    return(2.0*aa_*x + bb_); // TODO: should be (2*aa*x + bb) ?
  }
}

// d^2 LogStar(x)/dx^2
inline double el::InnerEL::LogStar2(double x) {
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
inline el::InnerEL::InnerEL(){}

/**
 * @brief Constructor for InnerEL with dimensions of G matrix as inputs for memory allocation.
 * 
 * @param n_obs    Number of observations.
 * @param n_eqs    Number of estimating equations.
 */
inline el::InnerEL::InnerEL(int n_obs, int n_eqs) {
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
  G_ = MatrixXd::Zero(n_eqs_,n_obs1_); // NEW: JAN 1
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
 * @param support    Whether to have support correction.
 * @param supp_a       Tuning parameter for support correction (referred as "a" in chen-et-al08).
 * @param lambda0    Initial value for lambda.
 */
inline void el::InnerEL::set_opts(const int& max_iter, const double& rel_tol, 
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
inline void el::InnerEL::set_opts(const int& max_iter, const double& rel_tol, 
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
 * @param support    Whether to have support correction.
 * @param supp_a       Tuning parameter for support correction (referred as "a" in chen-et-al08).
 */
inline void el::InnerEL::set_opts(const int& max_iter, const double& rel_tol, 
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
 * @param support    Whether to have support correction.
 */
inline void el::InnerEL::set_opts(const int& max_iter, const double& rel_tol, 
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
 * @param support    Whether to have support correction.
 * @param supp_a       Tuning parameter for support correction (referred as "a" in chen-et-al08).
 */
inline void el::InnerEL::set_opts(const bool& supp, const double& supp_a) {
  supp_ = supp;
  supp_a_ = supp_a;
  n_obs2_ = n_obs_+supp_;
}

/**
 * @brief Set support correction option only. 
 * 
 * @param support    Whether to have support correction.
 */
inline void el::InnerEL::set_opts(const bool& supp) {
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
inline void el::InnerEL::LambdaNR(int& n_iter, double& max_iter) {
  
  lambda_old_ = lambda0_; // set to initial value
  lambda_new_.fill(0.0); // may not be needed here..
  
  // Note: these two cannot be preallocate untill `support` is set
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
inline void el::InnerEL::EvalOmegas() {
  // G and lambdaNew must have been assigned
  if (lambda_new_ != lambda_new_) { // if lambdaNew is NaN 
    for (int ii=0; ii<n_obs2_; ii++) {
      omegas_(ii) = std::numeric_limits<double>::quiet_NaN();
    }
  }
  else {
    // if (supp_) adj_G(G_,supp_a_);
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
inline double el::InnerEL::LogEL() {
  // if omegas are NaN, return -Inf
  if (omegas_.head(n_obs2_) != omegas_.head(n_obs2_)) return -INFINITY;
  else return(omegas_.head(n_obs2_).array().log().sum());
}

// setters

/**
 * @brief Set the value of lambda (e.g. to be used directly to calculate omegas).
 */
inline void el::InnerEL::set_lambda(const Ref<const VectorXd>& lambda) {
  // lambda_old_ = lambda;
  lambda_new_ = lambda;
}

/**
 * @brief Set the value of omegas (e.g. to be used directly to calculate log EL).
 */
inline void el::InnerEL::set_omegas(const Ref<const VectorXd>& omegas) {
  omegas_.head(n_obs2_) = omegas; 
}

/**
 * @brief Set the value of G (e.g. to be used directly to calculate lambda or log EL).
 */
inline void el::InnerEL::set_G(const Ref<const MatrixXd>& G) {
  G_.block(0,0,n_eqs_,n_obs_) = G;
  if (supp_) adj_G(G_,supp_a_);
}

// getters

/**
 * @brief Get the value of lambda.
 */
inline VectorXd el::InnerEL::get_lambda() {
  return(lambda_new_);
}

/**
 * @brief Get the value of omegas.
 */
inline VectorXd el::InnerEL::get_omegas() {
    return(omegas_.head(n_obs2_));
}

/**
 * @brief Get the value of G.
 */
inline MatrixXd el::InnerEL::get_G() {
  return(G_.block(0,0,n_eqs_,n_obs2_));
}

/**
 * @brief Get the reference of G.
 */
inline Ref<MatrixXd> el::InnerEL::get_ref_G() {
  return Ref<MatrixXd>(G_.block(0,0,n_eqs_,n_obs2_));
}

/* TAKE OUT ALL MCMC SAMPLERS FOR NOW:

// This works for location models but also for multiple quantile levels
template<typename ELModel>
inline MatrixXd InnerEL<ELModel>::postSample(int nsamples, int nburn,
                MatrixXd ThetaInit, const Ref<const MatrixXd>& MwgSds,
                MatrixXd &RvDoMcmc, MatrixXd &Paccept) {
  int nThe = ThetaInit.rows(); // dimension of Theta
  int numThe = ThetaInit.cols(); // numer of Thetas
  MatrixXd Theta_chain(nThe*numThe,nsamples);
  MatrixXd ThetaOld = ThetaInit;
  MatrixXd ThetaNew = ThetaOld;
  MatrixXd ThetaProp = ThetaOld;
  ELModel::evalG(ThetaOld);
  int n_iter;
  double max_iter;
  LambdaNR(n_iter, max_iter);
  // TODO: what if not converged ?
  if (n_iter == max_iter && max_iter > rel_tol) {
    std::cout << "ThetaInit not valid." << std::endl;
    // return NULL;
  }
  EvalOmegas();
  double logELOld = LogEL(); 
  double logELProp;
  bool satisfy;
  double u;
  double a;
  double ratio;
  Paccept = MatrixXd::Zero(ThetaInit.rows(),ThetaInit.cols()); // TODO: need to initialize to 0 ?
  
  bool go_next;
  for (int ii=-nburn; ii<nsamples; ii++) {
    go_next = false;
    for (int kk=0; kk<numThe; kk++) {  
      if (go_next == true) break;
      for (int jj=0; jj<nThe; jj++) {
        if (RvDoMcmc(jj,kk)) {
          ThetaProp = ThetaOld;
          ThetaProp(jj,kk) += MwgSds(jj,kk)*R::norm_rand();
          // check if proposed Theta satisfies the constraint
          satisfy = false;
          ELModel::evalG(ThetaProp); // NEW: change G with ThetaProp
          LambdaNR(n_iter, max_iter);
          if (n_iter < max_iter) satisfy = true;
          // if does not satisfy, keep the old Theta
          if (satisfy == false) {
            go_next = true; // break out two loops
            break;
          }
          // if does satisfy
          u = R::unif_rand();
          // use the lambda calculate just now to get the logEL for Prop
          // to avoid an extra call of LambdaNR
          VectorXd logomegahat = 1/(1-(lambdaNew.transpose()*G).array());
          logomegahat = log(logomegahat.array()) - log(logomegahat.sum());
          // VectorXd logomegahat = log(1/(1-(lambdaNew.transpose()*ELModel::G).array())) -
          //   log((1/(1-(lambdaNew.transpose()*ELModel::G).array())).sum());
          logELProp = logomegahat.sum();
          ratio = exp(logELProp-logELOld);
          a = std::min(1.0,ratio);
          if (u < a) { // accepted
            Paccept(jj,kk) += 1; 
            ThetaNew = ThetaProp;
            ThetaOld = ThetaNew;
            logELOld = logELProp; // NEW: store the new one
          }
        }
      }
    }
    if (ii >= 0) {
      // Theta_chain.col(ii) = ThetaNew;
      Theta_chain.col(ii) = Map<VectorXd>(ThetaNew.data(), ThetaNew.size());
    }
  }
  Paccept /= (nsamples+nburn); 
  return(Theta_chain);
}

// mwgStep updates the idx entry of thetaCur
template<typename ELModel>
inline void InnerEL<ELModel>::mwgStep(VectorXd &thetaCur,
                                      const int &idx,
                                      const double &mwgsd,
                                      bool &accept, 
                                      double &logELCur) {
  int nThe = thetaCur.size();
  accept = false;
  VectorXd thetaProp = thetaCur;
  thetaProp(idx) += mwgsd*R::norm_rand();
  // sig2 has to be positive
  if (idx == nBet+nGam && thetaProp(idx) < 0) return;
  
  if (nThe == nBet) {
    // location mean regression and quantile regression (single quantile)
    ELModel::evalG(thetaProp);
  }
  else if (nThe == nBet + nGam + 1){
    // location-scale mean regression
    ELModel::evalG(thetaProp.head(nBet), 
                   thetaProp.segment(nBet,nGam), 
                   thetaProp.tail(1)(0),
                   VectorXd::Zero(0));
  }
  else {
    // location-scale quantile regression (single or multiple quantile)
    // using ELModel::nQts;
    ELModel::evalG(thetaProp.head(nBet), 
                   thetaProp.segment(nBet,nGam), 
                   thetaProp.segment(nBet+nGam,1)(0),
                   thetaProp.tail(nQts));
  }
  int n_iter;
  double max_iter;
  LambdaNR(n_iter, max_iter);
  bool satisfy = false;
  if (n_iter < max_iter || max_iter <= rel_tol) satisfy = true;
  // if does not satisfy, keep the old theta
  if (satisfy == false) return;
  // if does satisfy, might accept new theta
  double u = R::unif_rand();
  VectorXd logomegahat = log(1/(1-(lambdaNew.transpose()*ELModel::G).array())) -
    log((1/(1-(lambdaNew.transpose()*ELModel::G).array())).sum());
  double logELProp = logomegahat.sum();
  double ratio = exp(logELProp-logELCur);
  double a = std::min(1.0,ratio);
  if (u < a) { // accepted
    accept = true;
    thetaCur = thetaProp;
    logELCur = logELProp;
  }
}

template<typename ELModel>
inline MatrixXd InnerEL<ELModel>::postSampleAdapt(int nsamples, int nburn,
                                                  VectorXd thetaInit,
                                                  double *mwgSd, bool *rvDoMcmc,
                                                  VectorXd &paccept) {
  
  int nThe = thetaInit.size();
  MwgAdapt tuneMCMC(nThe, rvDoMcmc);
  bool *isAccepted = new bool[nThe];
  for (int ii=0; ii<nThe; ii++) {
    isAccepted[ii] = false;
  }
  MatrixXd theta_chain(nThe,nsamples);
  paccept = VectorXd::Zero(nThe);
  VectorXd thetaCur = thetaInit;
  if (nThe == nBet) {
    ELModel::evalG(thetaCur);
  }
  else if (nThe == nBet + nGam + 1){
    ELModel::evalG(thetaCur.head(nBet), 
                   thetaCur.segment(nBet,nGam), 
                   thetaCur.tail(1)(0),
                   VectorXd::Zero(0));
  }
  else {
    // using ELModel::nQts;
    ELModel::evalG(thetaCur.head(nBet), 
                   thetaCur.segment(nBet,nGam), 
                   thetaCur.segment(nBet+nGam,1)(0),
                   thetaCur.tail(nQts));
    // std::cout << "G = \n" << G << std::endl;
  }
  int n_iter;
  double max_iter;
  LambdaNR(n_iter, max_iter);
  // TODO: throw an error ??
  if (n_iter == max_iter && max_iter > rel_tol) {
    std::cout << "thetaInit not valid." << std::endl;
  }
  EvalOmegas();
  double logELCur = LogEL();
  // MCMC loop
  for(int ii=-nburn; ii<nsamples; ii++) {
    for(int jj=0; jj<nThe; jj++) {
      if(rvDoMcmc[jj]) {
        // std::cout << "isAccepted[" << jj << "] = " << isAccepted[jj] << std::endl;
        // modifies thetaCur's jj-th entry
        mwgStep(thetaCur,jj,mwgSd[jj],isAccepted[jj],logELCur);
        if (isAccepted[jj]) paccept(jj) += 1; // add 1 to paccept if accepted
      }
    }
    if (ii >= 0) {
      theta_chain.col(ii) = thetaCur;
    }
    tuneMCMC.adapt(mwgSd, isAccepted);
  }
  paccept /= (nsamples+nburn);
  delete[] isAccepted; // deallocate memory
  return(theta_chain);
}

*/

#endif
