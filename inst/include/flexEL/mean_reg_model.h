/**
 * @file mean_reg_model.h
 * 
 * @brief A base class for template classes InnerEL and InnerELC.
 */

#ifndef MEAN_REG_MODEL_H
#define MEAN_REG_MODEL_H

using namespace Eigen;

/* --------------------------------------------------------------------------- */

/**
* @brief EL namespace
* 
* Wrap the exported library components into a namespace called \b el to avoid potential naming conflicts with other libraries or user-defined headers.
*/
namespace flexEL {

/**
* @file       mean_reg_model.h
*
* @class      MeanRegModel
*
* @brief      A base class for template classes InnerEL and InnerELC to calculate estimating equations for mean regressions.
*/
class MeanRegModel {
  
  private:
    
    RowVectorXd yXb_; // y-X*beta
    RowVectorXd eZg_; // exp(-Z*gamma)
    RowVectorXd yXbeZg_; // (y-X*beta)*exp(-Z*gamma)
    RowVectorXd yXbeZg2_; // (y-X*beta)^2*exp(-2*Z*gamma)
    MatrixXd tG_; // t(G)
    int n_bet_, n_gam_; // dimensions of beta and gamma
    VectorXd y_; // vector of responses
    MatrixXd X_; // covariate matrix used in location function
    MatrixXd Z_; // covariate matrix used in scale function
    
  protected:
    
    int n_obs_; /**< number of observations (number of columns of G) */
    int n_eqs_; /**< number of estimating equations (number of rows of G) */
  
  public:
    
    // constructors
    MeanRegModel();
    MeanRegModel(int n_obs, int n_eqs);
    MeanRegModel(const Ref<const VectorXd>& y,
                 const Ref<const MatrixXd>& X);
    MeanRegModel(const Ref<const VectorXd>& y,
                 const Ref<const MatrixXd>& X,
                 const Ref<const MatrixXd>& Z);
    
    // set data functions
    void set_data(const Ref<const VectorXd>& y,
                  const Ref<const MatrixXd>& X);
    void set_data(const Ref<const VectorXd>& y,
                  const Ref<const MatrixXd>& X,
                  const Ref<const MatrixXd>& Z);
    
    // evaluate G matrix
    void eval_G(Ref<MatrixXd> G, 
                const Ref<const VectorXd>& beta);
    void eval_G(Ref<MatrixXd> G,
                const Ref<const VectorXd>& beta,
                const Ref<const VectorXd>& gamma,
                const double& sig2);
    
    // evaluate dGdt 
    // (p+q+1) x (p+q+1) x n_obs tensor in a (n_obs x (p+q+1)) x (p+q+1) matrix
    void eval_dGdt(Ref<MatrixXd> dGdt,
                   const Ref<const VectorXd>& beta);
    void eval_dGdt(Ref<MatrixXd> dGdt,
                   const Ref<const VectorXd>& beta,
                   const Ref<const VectorXd>& gamma,
                   const double& sig2);
    
    // get functions
    int get_n_obs();
    int get_n_eqs();
    
  };
}

/* --------------------------------------------------------------------------- */

// public functions

// ctors
/**
* @brief Default constructor for MeanRegModel.
*/
inline flexEL::MeanRegModel::MeanRegModel(){}

/**
* @brief Constructor for MeanRegModel with dimensions as inputs.
* 
* @param n_obs    Number of observations.
* @param n_eqs    Number of estimating equations.
*/
inline flexEL::MeanRegModel::MeanRegModel(int n_obs, int n_eqs) {
  n_obs_ = n_obs;
  n_eqs_ = n_eqs; // X gets passed as n_bet x n_obs matrix
  // TODO: pre-allocate space for data?
  // y_ = VectorXd::Zero(n_obs_);
}

/**
 * @brief Constructor for MeanRegModel with data as inputs (location model).
 * 
 * @param[in] y      Responses of length <code>n_obs</code>.
 * @param[in] X      Covariate matrix of dimension <code>n_bet</code> x <code>n_obs</code>.
 */
inline flexEL::MeanRegModel::MeanRegModel(const Ref<const VectorXd>& y,
                                          const Ref<const MatrixXd>& X) {
  y_ = y;
  X_ = X;
  n_obs_ = y.size();
  n_bet_ = X.rows(); // X gets passed as n_bet x n_obs matrix
  n_gam_ = 0;
  n_eqs_ = n_bet_;
}

/**
 * @brief Constructor for MeanRegModel with data as inputs (location-scale model).
 * 
 * @param[in] y      Responses of length <code>n_obs</code>.
 * @param[in] X      Covariate matrix of dimension <code>n_bet</code> x <code>n_obs</code>.
 * @param[in] Z      Covariate matrix of dimension <code>n_gam</code> x <code>n_obs</code>.
 */
inline flexEL::MeanRegModel::MeanRegModel(const Ref<const VectorXd>& y,
                                          const Ref<const MatrixXd>& X,
                                          const Ref<const MatrixXd>& Z) {
  y_ = y;
  X_ = X;
  Z_ = Z;
  n_obs_ = y.size();
  n_bet_ = X.rows(); // X gets passed as n_bet x n_obs matrix
  n_gam_ = Z.rows(); // Z gets passed as n_gam x n_obs matrix
  n_eqs_ = n_bet_+n_gam_+1;
}

// set_data location model (with default ctor)
/**
* @brief Set data for mean regression location model.
* 
* @param[in] y      Responses of length <code>n_obs</code>.
* @param[in] X      Covariate matrix of dimension <code>n_bet</code> x <code>n_obs</code>.
*/
inline void flexEL::MeanRegModel::set_data(const Ref<const VectorXd>& y,
                                           const Ref<const MatrixXd>& X) {
  y_ = y;
  X_ = X;
  n_obs_ = y.size();
  n_bet_ = X.rows(); // X gets passed as n_bet x n_obs matrix
  n_gam_ = 0;
  n_eqs_ = n_bet_;
}

// set_data location-scale model (with default ctor)
/**
* @brief Set data for mean regression location-scale model.
* 
* @param[in] y      Responses of length <code>n_obs</code>.
* @param[in] X      Covariate matrix of dimension <code>n_bet</code> x <code>n_obs</code>.
* @param[in] Z      Covariate matrix of dimension <code>n_gam</code> x <code>n_obs</code>.
*/
inline void flexEL::MeanRegModel::set_data(const Ref<const VectorXd>& y,
                                           const Ref<const MatrixXd>& X,
                                           const Ref<const MatrixXd>& Z) {
  y_ = y;
  X_ = X;
  Z_ = Z;
  n_obs_ = y.size();
  n_bet_ = X.rows(); // X gets passed as n_bet x n_obs matrix
  n_gam_ = Z.rows(); // Z gets passed as n_gam x n_obs matrix
  n_eqs_ = n_bet_+n_gam_+1;
}

// form the G matrix for location linear regression model
/**
* @brief Evaluate G matrix for mean regression location model.
* 
* @param[in] G        Matrix of dimension <code>n_eqs_ x n_obs_</code> where the calculated result is saved.
* @param[in] beta     Coefficient vector of length <code>n_bet_</code> in linear location function.
*/
inline void flexEL::MeanRegModel::eval_G(Ref<MatrixXd> G, const Ref<const VectorXd>& beta) {
  yXb_ = y_.transpose() - beta.transpose() * X_;
  tG_ = X_.transpose();
  tG_.array().colwise() *= yXb_.transpose().array();
  G = tG_.transpose();
}

// form the G matrix for location-scale linear regression model
/**
* @brief Evaluate G matrix for mean regression location-scale model.
*
* @param[in] G        Matrix of dimension <code>n_eqs_ x n_obs_</code> where the calculated result is saved.
* @param[in] beta     Coefficient vector of length <code>n_bet</code> in linear location function.
* @param[in] gamma    Coefficient vector of length <code>n_gam</code> in exponential scale function.
* @param[in] sig2     Scale parameter in scale function.
*/
inline void flexEL::MeanRegModel::eval_G(Ref<MatrixXd> G, 
                                         const Ref<const VectorXd>& beta,
                                         const Ref<const VectorXd>& gamma,
                                         const double &sig2) {
  eZg_.array() = (-gamma.transpose()*Z_).array().exp();
  yXbeZg_.array() = (y_.transpose()-beta.transpose()*X_).array() * eZg_.array();
  yXbeZg2_.array() = yXbeZg_.array()*yXbeZg_.array();
  
  tG_ = MatrixXd::Zero(n_obs_, n_eqs_); // NEW: DEC 25
  tG_.block(0,0,n_obs_,n_bet_) = X_.transpose();
  tG_.block(0,0,n_obs_,n_bet_).array().colwise() *= yXbeZg_.transpose().array() * eZg_.transpose().array();
  tG_.block(0,n_bet_,n_obs_,n_gam_) = Z_.transpose();
  tG_.block(0,n_bet_,n_obs_,n_gam_).array().colwise() *= (1.0-yXbeZg2_.transpose().array());
  tG_.rightCols(1).array() = 1/sig2*yXbeZg2_.transpose().array()-1;
  G = tG_.transpose();
}

/**
 * @brief Evaluate derivative of G matrix w.r.t beta for mean regression location model.
 * 
 * @param[in] G        Matrix of dimension <code>n_eqs_ x n_obs_</code> where the calculated result is saved.
 * @param[in] beta     Coefficient vector of length <code>n_bet_</code> in linear location function.
 */
inline void flexEL::MeanRegModel::eval_dGdt(Ref<MatrixXd> dGdt, 
                                            const Ref<const VectorXd>& beta) {
  for(int ii=0; ii<n_obs_; ii++) {
    dGdt.block(ii*n_bet_,0,n_bet_,n_bet_) = -X_.col(ii)*X_.col(ii).transpose();
  }
}

/**
 * @brief Evaluate derivative of G matrix w.r.t theta for mean regression location-scale model.
 *
 * @param[in] G        Matrix of dimension <code>n_eqs_ x n_obs_</code> where the calculated result is saved.
 * @param[in] beta     Coefficient vector of length <code>n_bet</code> in linear location function.
 * @param[in] gamma    Coefficient vector of length <code>n_gam</code> in exponential scale function.
 * @param[in] sig2     Scale parameter in scale function.
 */
inline void flexEL::MeanRegModel::eval_dGdt(Ref<MatrixXd> dGdt,
                                            const Ref<const VectorXd>& beta,
                                            const Ref<const VectorXd>& gamma,
                                            const double& sig2) {
  eZg_.array() = (-gamma.transpose()*Z_).array().exp();
  yXb_ = y_.transpose() - beta.transpose() * X_;
  // std::cout << "eZg_ = " << eZg_.transpose() << std::endl;
  // std::cout << "yXb_ = " << yXb_.transpose() << std::endl;
  // std::cout << "yXb_*yXb_ = " << yXb_*yXb_ << std::endl;
  // std::cout << "sig2*sig2 = " << sig2*sig2 << std::endl;
  double eZg2_ii;
  double yXb_ii;
  for(int ii=0; ii<n_obs_; ii++) {
    int br_start = ii*(n_bet_+n_gam_+1);
    eZg2_ii = eZg_(ii)*eZg_(ii);
    yXb_ii = yXb_(ii);
    // std::cout << "eZg2_ii = " << eZg2_ii << std::endl;
    // std::cout << "yXb_ii = " << yXb_ii << std::endl;
    // dGi/db
    dGdt.block(br_start,0,n_bet_,n_bet_) = -eZg2_ii*X_.col(ii)*X_.col(ii).transpose();
    dGdt.block(br_start+n_bet_,0,n_gam_,n_bet_) = 2/sig2*eZg2_ii*yXb_ii*Z_.col(ii)*X_.col(ii).transpose();
    dGdt.block(br_start+n_bet_+n_gam_,0,1,n_bet_) = -2/sig2*eZg2_ii*yXb_ii*X_.col(ii).transpose();
    // dGi/dg
    dGdt.block(br_start,n_bet_,n_bet_,n_gam_) = -2*eZg2_ii*yXb_ii*X_.col(ii)*Z_.col(ii).transpose();
    dGdt.block(br_start+n_bet_,n_bet_,n_gam_,n_gam_) = 2/sig2*eZg2_ii*yXb_ii*yXb_ii*Z_.col(ii)*Z_.col(ii).transpose();
    dGdt.block(br_start+n_bet_+n_gam_,n_bet_,1,n_gam_) = -2/sig2*eZg2_ii*yXb_ii*yXb_ii*Z_.col(ii).transpose();
    // dGi/ds2
    dGdt.block(br_start,n_bet_+n_gam_,n_bet_,1) = VectorXd::Zero(n_bet_);
    dGdt.block(br_start+n_bet_,n_bet_+n_gam_,n_gam_,1) = 1/(sig2*sig2)*eZg2_ii*yXb_ii*yXb_ii*Z_.col(ii);
    // std::cout << "dGi/ds2 = " << -1/(sig2*sig2)*eZg2_ii*yXb_ii*yXb_ii << std::endl;
    // std::cout << "dGdt.block(br_start+n_bet_+n_gam_,n_bet_+n_gam_,1,1) = " << dGdt.block(br_start+n_bet_+n_gam_,n_bet_+n_gam_,1,1) << std::endl;
    dGdt.block(br_start+n_bet_+n_gam_,n_bet_+n_gam_,1,1).array() = -1/(sig2*sig2)*eZg2_ii*yXb_ii*yXb_ii;
  }
}

/**
 * @brief Return the number of observations.
 */
inline int flexEL::MeanRegModel::get_n_obs() {
  return n_obs_;
}

/**
 * @brief Return the number of estimating equations.
 */
inline int flexEL::MeanRegModel::get_n_eqs() {
  return n_eqs_;
}

#endif
