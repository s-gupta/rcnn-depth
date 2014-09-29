/*
 * Functors for common distance metrics on matrices.
 */
#ifndef MATH__MATRICES__FUNCTORS__MATRIX_DISTANCE_FUNCTORS_HH
#define MATH__MATRICES__FUNCTORS__MATRIX_DISTANCE_FUNCTORS_HH

#include "functors/distanceable_functors.hh"
#include "math/math.hh"
#include "math/matrices/matrix.hh"

namespace math {
namespace matrices {
namespace functors {
/*
 * Imports.
 */
using ::functors::distanceable_functor;
using math::matrices::matrix;
using math::abs;
using math::pow;
using math::sqrt;

/*
 * L1 distance between matrices.
 */
template <typename T = double>
class matrix_L1_distance : public distanceable_functor<const matrix<T>, T> {
public:
   T operator()(const matrix<T>& m0, const matrix<T>& m1) const {
      matrix<T>::assert_dims_equal(m0._dims, m1._dims);
      T dist = T();
      for (unsigned long n = 0; n < m0._size; n++)
         dist += abs(m1._data[n] - m0._data[n]);
      return dist;
   }
};

/*
 * L2 distance between matrices.
 */
template <typename T = double>
class matrix_L2_distance : public distanceable_functor<const matrix<T>, T> {
public:
   T operator()(const matrix<T>& m0, const matrix<T>& m1) const {
      matrix<T>::assert_dims_equal(m0._dims, m1._dims);
      T dist = T();
      for (unsigned long n = 0; n < m0._size; n++) {
         T diff = m1._data[n] - m0._data[n];
         dist += diff*diff;
      }
      return sqrt(dist);
   }
};

/*
 * Ln distance between matrices.
 */
template <typename T = double>
class matrix_Ln_distance : public distanceable_functor<const matrix<T>, T> {
public:
   /*
    * Constructor.
    * Specify the value of n.
    */
   explicit matrix_Ln_distance(const T& n)
    : _NN(n)
   { }

   /*
    * Copy constructor.
    */
   explicit matrix_Ln_distance(const matrix_Ln_distance<T>& f)
    : _NN(f._NN)
   { }
   
   /*
    * Compute the Ln distance.
    */
   T operator()(const matrix<T>& m0, const matrix<T>& m1) const {
      matrix<T>::assert_dims_equal(m0._dims, m1._dims);
      T dist = T();
      for (unsigned long n = 0; n < m0._size; n++)
         dist += pow(m1._data[n] - m0._data[n], _NN);
      return pow(dist, T(1)/_NN);
   }
   
protected:
   T _NN;
};

/*
 * Chi-squared distance between matrices.
 */
template <typename T = double>
class matrix_X2_distance : public distanceable_functor<const matrix<T>, T> {
public:
   T operator()(const matrix<T>& m0, const matrix<T>& m1) const {
      matrix<T>::assert_dims_equal(m0._dims, m1._dims);
      T dist = T();
      for (unsigned long n = 0; n < m0._size; n++) {
         T diff = m1._data[n] - m0._data[n];
         T sum  = m1._data[n] + m0._data[n];
         if (diff != T())
            dist += diff*diff / sum;
      }
      return dist/T(2);
   }
};

/*
 * Globally accessible set of useful matrix distance functors.
 */
template <typename T = double>
class matrix_distance_functors {
public:
   static const matrix_L1_distance<T>& L1_distance();
   static const matrix_L2_distance<T>& L2_distance();
   static const matrix_X2_distance<T>& X2_distance();
};

template <typename T>
const matrix_L1_distance<T>& matrix_distance_functors<T>::L1_distance() {
   static const matrix_L1_distance<T>* f = new matrix_L1_distance<T>();
   return *f;
}

template <typename T>
const matrix_L2_distance<T>& matrix_distance_functors<T>::L2_distance() {
   static const matrix_L2_distance<T>* f = new matrix_L2_distance<T>();
   return *f;
}

template <typename T>
const matrix_X2_distance<T>& matrix_distance_functors<T>::X2_distance() {
   static const matrix_X2_distance<T>* f = new matrix_X2_distance<T>();
   return *f;
}

} /* namespace functors */
} /* namespace matrices */
} /* namespace math */

#endif
