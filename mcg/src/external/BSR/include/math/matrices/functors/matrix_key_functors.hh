/*
 * Functors for comparing matrices to axis-aligned half-space keys.
 */
#ifndef MATH__MATRICES__FUNCTORS__MATRIX_KEY_FUNCTORS_HH
#define MATH__MATRICES__FUNCTORS__MATRIX_KEY_FUNCTORS_HH

#include "collections/abstract/collection.hh"
#include "interfaces/kd_treeable.hh"
#include "functors/kd_treeable_functors.hh"
#include "functors/keyable_functors.hh"
#include "lang/exceptions/ex_invalid_argument.hh"
#include "lang/iterators/iterator.hh"
#include "lang/pointers/auto_ptr.hh"
#include "math/matrices/matrix.hh"

namespace math {
namespace matrices {
namespace functors {
/*
 * Imports.
 */
using collections::abstract::collection;
using interfaces::kd_tree_key;
using ::functors::kd_treeable_key_distance_functor;
using ::functors::keyable_compare_functor;
using ::functors::keyable_split_functor;
using lang::exceptions::ex_invalid_argument;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;
using math::matrices::matrix;

/*
 * Compare matrix element at dimension specified in key to key value.
 */
template <typename T = double>
class matrix_key_compare
 : public keyable_compare_functor< const matrix<T>, const kd_tree_key<T> > {
public:
   int operator()(const matrix<T>& m, const kd_tree_key<T>& k) const {
      unsigned long split_dim = k.split_dimension();
      T split_val = k.split_value();
      T val = m[split_dim];
      return ((val < split_val) ? (-1) : ((val > split_val) ? 1 : 0));
   }
};

/*
 * Compute key for splitting a collection of matrices (regarded as 
 * vectors in a possibly high dimensional feature space) along the 
 * dimension of the feature vector with highest variance.  Throw an 
 * invalid argument exception if the matrices do not all have the 
 * same number of elements or the collection is empty.
 */
template <typename T = double>
class matrix_key_split
 : public keyable_split_functor< matrix<T>, kd_tree_key<T> > {
public:
   kd_tree_key<T> operator()(const collection< matrix<T> >& c) const {
      /* check that collection is nonempty */
      unsigned long n_matrices = c.size();
      if (n_matrices == 0)
         throw ex_invalid_argument(
            "attempt to split empty collection"
         );
      /* check that all matrices have the same size */
      auto_ptr< iterator< matrix<T> > > i = c.iter_create();
      unsigned long n_dims = i->next().size();
      while (i->has_next()) {
         if (i->next().size() != n_dims)
            throw ex_invalid_argument(
               "all matrices in collection must have the same size"
            );
      }
      i.reset();
      /* compute variance along each dimension */
      matrix<T> sum_sq_x(n_dims, 1);
      matrix<T> sum_x_sq(n_dims, 1);
      matrix<T> variance(n_dims, 1);
      i = c.iter_create();
      while (i->has_next()) {
         matrix<T>& m = i->next();
         for (unsigned long n = 0; n < n_dims; n++) {
            T val = m[n];
            sum_sq_x[n] += val;
            sum_x_sq[n] += val*val;
         }
      }
      i.reset();
      for (unsigned long n = 0; n < n_dims; n++) {
         sum_x_sq[n] /= T(n_matrices);
         sum_sq_x[n] /= T(n_matrices);
         sum_sq_x[n] *= sum_sq_x[n];
         variance[n] = sum_x_sq[n] - sum_sq_x[n];
      }
      /* find dimension with maximum variance */
      unsigned long max_var_dim = 0;
      T max_var = variance[0];
      for (unsigned long n = 0; n < n_dims; n++) {
         if (variance[n] > max_var) {
            max_var = variance[n];
            max_var_dim = n;
         }
      }
      /* extract coordinates along this dimension */
      matrix<T> vals(n_matrices, 1);
      i = c.iter_create();
      for (unsigned long n = 0; n < n_matrices; n++) {
         matrix<T>& m = i->next();
         vals[n] = m[max_var_dim];
      }
      i.reset();
      /* find median of coordinates in this dimension */
      vals.sort();
      T min_val   = vals[0];
      T split_val = vals[(n_matrices/2)];
      T max_val   = vals[n_matrices-1];
      /* make sure at least one element falls on left side */
      if (min_val == split_val)
         split_val = max_val;
      return kd_tree_key<T>(max_var_dim, split_val);
   }
};

/*
 * Compute distance from matrix (viewed as a vector) to
 * the axis-aligned half-space represented by the key.
 */
template <typename T = double>
class matrix_key_distance
 : public kd_treeable_key_distance_functor<const matrix<T>,T> {
public:
   T operator()(const matrix<T>& m, const kd_tree_key<T>& k) const {
      T dist = k.split_value() - m[k.split_dimension()];
      return ((dist < T()) ? (-dist) : dist);
   }
};
 
/*
 * Globally accessible set of matrix key functors.
 */
template <typename T = double>
class matrix_key_functors {
public:
   static const matrix_key_compare<T>&  f_key_compare();
   static const matrix_key_split<T>&    f_key_split();
   static const matrix_key_distance<T>& f_key_distance();
};

template <typename T>
const matrix_key_compare<T>& matrix_key_functors<T>::f_key_compare() {
   static const matrix_key_compare<T>* f = new matrix_key_compare<T>();
   return *f;
}

template <typename T>
const matrix_key_split<T>& matrix_key_functors<T>::f_key_split() {
   static const matrix_key_split<T>* f = new matrix_key_split<T>();
   return *f;
}

template <typename T>
const matrix_key_distance<T>& matrix_key_functors<T>::f_key_distance() {
   static const matrix_key_distance<T>* f = new matrix_key_distance<T>();
   return *f;
}

} /* namespace functors */
} /* namespace matrices */
} /* namespace math */

#endif
