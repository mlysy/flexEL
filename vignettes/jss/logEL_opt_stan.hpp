/// @file logEL_opt_stan.hpp
/// @brief Example of wrapping the `logEL()` and its gradient in Stan.

#ifndef stanEL_logEL_opt_stan_hpp
#define stanEL_logEL_opt_stan_hpp 1

// this line tells the compiler to look for files in the installed flexEL
//[[Rcpp::depends(flexEL)]]
// this lines tells the compiler specifically which file to look for
#include <flexEL/gen_el.h>
// link to the stan math library
#include <stan/math.hpp>
#include <Rcpp.h>

template <typename T_G, typename T_m, typename T_r,
          typename T_s, typename T_d>
inline auto logEL_opt(const Eigen::Matrix<T_G,
		      Eigen::Dynamic,
		      Eigen::Dynamic>& G,
		      const T_m& max_iter,
		      const T_r& rel_tol,
		      const T_s& supp_adj,
		      const T_d& dnc_inf) {

  // import stan names
  using stan::math::operands_and_partials;
  using stan::partials_return_t;
  // create the return type
  using T_partials_return = partials_return_t<T_G, T_m, T_r, T_s, T_d>;
  // create the gradient wrt Eigen::Matrix type (i.e., stan matrix)
  using T_partials_matrix = Eigen::Matrix<T_partials_return,
                                          Eigen::Dynamic, Eigen::Dynamic>;

  // Dimensions of input G matrix
  //  input for flexEL::InnerEL
  int n_obs = G.cols();
  int n_eqs = G.rows();

  // create the return object containing the function and its gradients
  operands_and_partials<Eigen::Matrix<T_G, Eigen::Dynamic,
                                      Eigen::Dynamic>, T_r>
    ops_partials(G, rel_tol);
  // Evaluate the function and its gradient,
  // we'll only evaluate the gradient if necessary.
  bool return_dldG = !stan::is_constant_all<T_G>::value;
  T_partials_matrix dldG;
  // dummy variable for nonexistent gradient through rel_tol.
  T_partials_return dummy_zero = 0;
  // only allocate memory if necessary
  if(return_dldG) dldG = Eigen::MatrixXd::Zero(n_eqs, n_obs);
  double logel;

  // The function object:
  flexEL::GenEL IL(n_obs, n_eqs);
  //Specified values for options:
  int max_iter_ = value_of(max_iter);
  double rel_tol_ = value_of(rel_tol);
  bool supp_adj_ = value_of(supp_adj) > 0;
  bool dnc_inf_ = value_of(dnc_inf) > 0;
  IL.set_max_iter(max_iter_);
  IL.set_rel_tol(rel_tol_);
  IL.set_supp_adj(supp_adj_);

  // wrap input as Eigen::MatrixXd
  const auto& G_ = value_of(G);

  // Output:
  logel = IL.logel_grad(dldG, G_);
  if(dnc_inf_) {
    // return -Inf if NR has not converged
    int nr_iter;
    double nr_err;
    IL.get_diag(nr_iter, nr_err);
    if((nr_iter >= max_iter_) && (nr_err > rel_tol_)) {
      if(return_dldG) {
	dldG.setZero();
	ops_partials.edge1_.partials_ = dldG;
      }
      return ops_partials.build(negative_infinity());
    }
  }  
  if(return_dldG) {
    // put the gradient into the return object
    ops_partials.edge1_.partials_ = dldG;
  }

  if(!stan::is_constant_all<T_r>::value) {
    ops_partials.edge2_.partials_[0] = dummy_zero;
  }

  // put the function value into the return object
  T_partials_return ans(logel);
  return ops_partials.build(ans);
}

template <typename T0__, typename T2__>
typename boost::math::tools::promote_args<T0__, T2__>::type
logEL_opt(const Eigen::Matrix<T0__, Eigen::Dynamic, Eigen::Dynamic>& G,
	  const int& max_iter,
	  const T2__& rel_tol,
	  const int& supp_adj,
	  const int& dnc_inf, std::ostream* pstream__) {
  return logEL_opt(G, max_iter, rel_tol, supp_adj, dnc_inf);
}

#endif
