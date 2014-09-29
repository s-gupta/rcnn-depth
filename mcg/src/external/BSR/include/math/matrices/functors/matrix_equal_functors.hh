/*
 * Functors for checking equality of matrices.
 */
#ifndef MATH__MATRICES__FUNCTORS__MATRIX_EQUAL_FUNCTORS_HH
#define MATH__MATRICES__FUNCTORS__MATRIX_EQUAL_FUNCTORS_HH

#include "functors/equalable_functors.hh"
#include "math/matrices/matrix.hh"

namespace math {
namespace matrices {
namespace functors {
/*
 * Imports.
 */
using ::functors::equalable_functor;
using math::matrices::matrix;

/*
 * Exact equality of all elements.
 */
template <typename T = double>
class matrix_equal : public equalable_functor< const matrix<T> > {
public:
   bool operator()(const matrix<T>& m0, const matrix<T>& m1) const {
      matrix<T>::assert_dims_equal(m0._dims, m1._dims);
      bool eq = true;
      for (unsigned long n = 0; (n < m0._size) && eq; n++)
         eq = (m0._data[n] == m1._data[n]);
      return eq;
   }
};

/*
 * Globally accessible set of useful matrix eqaulity test functors.
 */
template <typename T = double>
class matrix_equal_functors {
public:
   static const matrix_equal<T>& f_equal();
};

template <typename T>
const matrix_equal<T>& matrix_equal_functors<T>::f_equal() {
   static const matrix_equal<T>* f = new matrix_equal<T>();
   return *f;
}

} /* namespace functors */
} /* namespace matrices */
} /* namespace math */

#endif
