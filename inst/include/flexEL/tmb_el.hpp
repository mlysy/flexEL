/// @file tmb_el.hpp
///
/// @brief TMB wrapper for `flexEL::GenEL<Type>.logel()`.
///
/**
    While `flexEL` is fully templated, we wouldn't want to use autodiff on the Newton-Raphson solver directly.  Rather, we would use `flexEL::GenEL<Type>.logel_grad()`.  Here's how this would work:
    
    1. `TMB_ATOMIC_VECTOR_FUNCTION` will create a function `logel_atomic()` with a custom reverse-mode rule for us.

    - The forward pass (`ATOMIC_DOUBLE`) has input and output of type `CppAD::vector<double>`:

        - The input (`tx`) should be a flattened version of `G`, prepended with the row/col dimensions and any other information, e.g., support correction, number of NR iterations, etc.

	- The output (`ty`) would normally be the result of `flexEL::GenEL<double>.logel(G)`, which is a scalar.  However, to efficiently compute the reverse-mode algorithm we actually need the intermediate results obtainted by calling `flexEL::GenEL<double>.logel_full()`, which returns `logel` as scalar and `omega` and `lambda` as intermediate vectors.  These vectors then need to be concatenated into a vector.  Thus, the `OUTPUT_DIM` macro needs to be set to the sum of the elements in `G`, `omega`, and `lambda`.

    - The backward pass (`ATOMIC_REVERSE`) consists of the body of a function with inputs and outputs all being of type `CppAD::vector<Type>`.

        -  The inputs are `tx`, `ty`, and `py`.  The first two correspond to the same quantities as in `ATOMIC_DOUBLE`, but having been converted to type `CppAD::vector<Type>`.  `py` is a vector of the same length as `ty`.

        - The output is `px`, a vector of the same length as `tx`.

	- The body of the backward pass first unpacks `G` from `tx`, and `omega` and `lambda` from `ty`.  Then it uses `flexEL::GenEL<Type>.logel_grad()` to compute `dldG`, the gradient of `logel_atomic()` with respect to `G`.  Finally, we set 
            ```
            px[i+n_nondiff] = d/dG[i] sum(logel(G) * py[0]) =  
                  = dldG[i] * py[0],
            ```
	    where `n_nondiff` is the number of nondifferentiable arguments prepended to the flattened version of `G` (i.e., dimensions of `G` and NR control parameters).

    2.  `logel_atomic()` is a function which takes vector inputs and outputs, each of which contain the flattened version of a matrix (`G` or `dldG`) along with some elements (corresponding to the nondifferentiable inputs) which will never get used.  So the final `logel()` function is just a wrapper which does the arguent packing and unpacking for us.

*/

#ifndef FLEXEL_TMB_EL_HPP
#define FLEXEL_TMB_EL_HPP

#include "gen_el.h"
// #include <TMB.hpp> // don't include since it's not header-guarded
#include <Eigen/Dense>
#include <Eigen/Core>

namespace flexEL {

  using namespace Eigen;

  // typedefs
  template <class Type>
  using cRefMatrix_t = const Ref <const Matrix<Type, Dynamic, Dynamic> >;
  template <class Type>
  using MatrixX_t = Matrix<Type, Dynamic, Dynamic>;
  template <class Type>
  using VectorX_t = Matrix<Type, Dynamic, 1>;

  // Sumber of nondifferentiable inputs.  See `tx_to_G()`.
  const int N_NONDIFF = 5;

  /// Copy contiguous elements of a `CppAD::vector<Type>` to a `VectorX_t<Type>`.
  ///
  /// The intented usage is something like
  /// ```
  /// const VectorX_t<Type> y = copy_vector<Type>(x)
  /// ```
  /// And this _does_ make a copy of the data: see <https://stackoverflow.com/questions/52261389/how-to-convert-an-stdvector-to-a-matrix-in-eigen>.
  ///
  /// @param[in] x A `CppAD::vector<Type>` object.
  /// @param[in] start The starting index.
  /// @param[in] size The number of elements to map.
  /// @return A `VectorX_t<Type>` consisting of the elements `x[start], ..., x[start + size-1]`.
  template <class Type>
  const VectorX_t<Type> copy_vector(const CppAD::vector<Type>& x,
				    int start, int size) {
    return Map<const VectorX_t<Type>,Unaligned>(x.data()+start, size);
  }

  /// Copy contiguous elements of a `CppAD::vector<Type>` to a `MatrixX_t<Type>`.
  ///
  /// @param[in] x A `CppAD::vector<Type>` object.
  /// @param[in] start The starting index.
  /// @param[in] rows The number of rows in the output matrix.
  /// @param[in] cols The number of columns in the output matrix.
  /// @return A `Matrix_t<Type>` object with dimensions `rows x cols` and consisting of the elements `x[start], ..., x[start + (row * cols)-1]`.
  template <class Type>
  const MatrixX_t<Type> copy_matrix(const CppAD::vector<Type>& x,
				    int start, int rows, int cols) {
    return Map<const MatrixX_t<Type>,Unaligned>(x.data()+start, rows, cols);
  }

  /// Extract G matrix from contiguous inputs.
  ///
  /// The contiguous input is assumed to start with the following scalars prior to the flattened `G`:
  ///
  /// - `n_eqs`, `n_obs`: The number of rows and columns of `G`.
  /// - `supp_adj`: Whether or not to do support adjustment.
  /// - `max_iter`: Maximum number of NR iterations.
  /// - `rel_tol`: Relative tolerance for NR convergence.
  ///
  /// @param[in] tx Vector of contiguous inputs.
  /// @return The matrix `G`.
  template <class Type>
  const MatrixX_t<Type> tx_to_G(const CppAD::vector<Type>& tx) {
    int n_eqs = CppAD::Integer(tx[0]);
    int n_obs = CppAD::Integer(tx[1]);
    return copy_matrix<Type>(tx, N_NONDIFF, n_eqs, n_obs);
  }

  /// Vectorized wrapper to `flexEL::GenEL<Type>.logel_full()`.
  ///
  /// Both inputs and outputs are contiguous vectors.  For the former, see `tx_to_G()`.  For the latter, the contiguous elements are:
  /// - `log_el`: The EL loglikelihood itself (scalar).
  /// - `has_conv`: Whether or not the NR converged (scalar).
  /// - `omega`, `lambda`: Auxillary outputs of `GenEL::logel_full()`.
  ///
  /// @param[out] ty Vector of continuous outputs.
  /// @param[in] tx Vector of contiguous inputs.
  template <class Type>
  void logel_vector(CppAD::vector<Type>& ty, const CppAD::vector<Type>& tx) {
    // unpack inputs from tx
    int n_eqs = int(tx[0]);
    int n_obs = int(tx[1]);
    bool supp_adj = bool(tx[2]);
    int n_obs2 = n_obs + supp_adj;
    int max_iter = int(tx[3]);
    Type rel_tol = Type(tx[4]);
    MatrixX_t<Type> G = tx_to_G(tx);
    // weights input needed for logel_full()
    VectorX_t<Type> weights = VectorX_t<Type>::Ones(n_obs);
    // outputs
    VectorX_t<Type> lambda(n_eqs);
    VectorX_t<Type> omega(n_obs2);
    // initialize GenEL object
    flexEL::GenEL<Type> gel(n_obs, n_eqs);
    gel.set_max_iter(max_iter);
    gel.set_rel_tol(rel_tol);
    gel.set_supp_adj(supp_adj);
    // calculations
    Type log_el = gel.logel_full(omega, lambda, G, weights);
    bool has_conv = gel.has_converged_nr();
    if(!has_conv) {
      log_el = -std::numeric_limits<double>::infinity();
    }
    // pack outputs into ty
    ty[0] = log_el;
    ty[1] = Type(has_conv);
    // NOTE: this makes a copy of the Eigen objects,
    // see <https://stackoverflow.com/questions/26094379/typecasting-eigenvectorxd-to-stdvector>
    VectorX_t<Type>::Map(&ty[2], n_obs2) = omega;
    VectorX_t<Type>::Map(&ty[2+n_obs2], n_eqs) = lambda;
    return;
  }

  TMB_ATOMIC_VECTOR_FUNCTION(
			     // ATOMIC_NAME
			     logel_atomic
			     ,
			     // OUTPUT_DIM
			     2 + CppAD::Integer(tx[0] + tx[1] + tx[2]),
			     // ATOMIC_DOUBLE
			     logel_vector<double>(ty, tx),
			     // ATOMIC_REVERSE
			     // unpack inputs from tx
			     int n_eqs = CppAD::Integer(tx[0]);
			     int n_obs = CppAD::Integer(tx[1]);
			     bool supp_adj = bool(CppAD::Integer(tx[2]));
			     int n_obs2 = n_obs + supp_adj;
			     int max_iter = CppAD::Integer(tx[3]);
			     Type rel_tol = Type(tx[4]);
			     MatrixX_t<Type> G = tx_to_G(tx);
			     // calculate gradient only if NR has converged
			     MatrixX_t<Type> dldG(n_eqs, n_obs);
			     bool has_conv = bool(CppAD::Integer(ty[1]));
			     if(!has_conv) {
			       dldG.setConstant(std::numeric_limits<Type>::quiet_NaN());
			     } else {
			       // unpack auxillary outputs from ty
			       int start = 2; // first two arguments are log_el and has_conv
			       VectorX_t<Type> omega = copy_vector<Type>(ty, start, n_obs2);
			       start += n_obs2;
			       VectorX_t<Type> lambda = copy_vector<Type>(ty, start, n_eqs);
			       // initialize GenEL object
			       flexEL::GenEL<Type> gel(n_obs, n_eqs);
			       gel.set_max_iter(max_iter);
			       gel.set_rel_tol(rel_tol);
			       gel.set_supp_adj(supp_adj);
			       // calculate gradient
			       Type sum_weights = Type(n_obs2);
			       gel.logel_grad(dldG, omega, lambda, sum_weights);
			     }
			     // compute the reverse-mode rule
			     for(int ii=0; ii<N_NONDIFF; ii++) {
			       px[ii] = Type(0.);
			     }
			     for(int jj=0; jj<n_obs; jj++) {
			       for(int ii=0; ii<n_eqs; ii++) {
				 px[jj*n_eqs + ii + N_NONDIFF] = dldG(ii,jj) * py[0];
			       }
			     }
			     )

  /// Wrapper to logel_atomic() for Eigen types.
  ///
  /// @param G Matrix of size `n_eqs x n_obs`.
  /// @param max_iter Maximum number of NR iterations.
  /// @param rel_tol Relative tolerance for NR convergence.
  /// @param supp_adj Whether or not to perform support adjustment.
  ///
  /// @return The return value of `flexEL::GenEL<Type>.logel()`.
			     template<class Type>
			     Type logel(cRefMatrix_t<Type>& G,
					int max_iter = 100,
					Type rel_tol = Type(1e-7),
					bool supp_adj = false) {
    CppAD::vector<Type> tx(G.size() + N_NONDIFF);
    int n_eqs = G.rows();
    int n_obs = G.cols();
    tx[0] = Type(n_eqs);
    tx[1] = Type(n_obs);
    tx[2] = Type(supp_adj);
    tx[3] = Type(max_iter);
    tx[4] = Type(rel_tol);
    for(int jj=0; jj<n_obs; jj++) {
      for(int ii=0; ii<n_eqs; ii++) {
	tx[jj*n_eqs + ii + N_NONDIFF] = G(ii,jj);
      }
    }
    return logel_atomic(tx)[0];
  }

} // end namespace flexEL

#endif
