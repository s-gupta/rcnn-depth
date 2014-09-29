/*
 * Matrix (thread-safe).
 */
#ifndef MATH__MATRICES__MATRIX_HH
#define MATH__MATRICES__MATRIX_HH

#include "concurrent/threads/synchronization/locks/auto_read_lock.hh"
#include "concurrent/threads/synchronization/locks/auto_read_read_lock.hh"
#include "concurrent/threads/synchronization/locks/auto_read_write_lock.hh"
#include "concurrent/threads/synchronization/locks/auto_write_lock.hh"
#include "concurrent/threads/synchronization/synchronizables/synchronizable.hh"
#include "concurrent/threads/synchronization/synchronizables/unsynchronized.hh"
#include "config/safety.hh"
#include "functors/comparable_functors.hh"
#include "io/serialization/serial_input_stream.hh"
#include "io/serialization/serial_output_stream.hh"
#include "io/serialization/serializers.hh"
#include "io/streams/ios.hh"
#include "io/streams/iomanip.hh"
#include "io/streams/ostream.hh"
#include "io/streams/ostringstream.hh"
#include "lang/array.hh"
#include "lang/exceptions/ex_index_out_of_bounds.hh"
#include "lang/exceptions/ex_invalid_argument.hh"
#include "math/complex.hh"
#include "math/math.hh"
#include "math/matrices/exceptions/ex_matrix_singular.hh"
#include "math/matrices/exceptions/ex_matrix_dimension_mismatch.hh"
#include "math/random/generators/rand_gen.hh"
#include "math/random/generators/rand_gen_uniform.hh"

/* 
 * Enable/disable bounds checking for matrices.
 */
#if CONFIG__SAFETY__CHECK_BOUNDS
   #define MATH__MATRICES__MATRIX__CHECK_BOUNDS (true)
#else
   #define MATH__MATRICES__MATRIX__CHECK_BOUNDS (false)
#endif

/*
 * Declare existence of libraries that will need friend access to matrices.
 */
namespace math {
namespace libraries {
   class lib_image;
   class lib_matrix;
   class lib_signal;
} /* namespace libraries */
} /* namespace math */

/*
 * Declare existence of functors that will need friend access to matrices.
 */
namespace math {
namespace matrices {
namespace functors {
   template <typename T> class matrix_L1_distance;
   template <typename T> class matrix_L2_distance;
   template <typename T> class matrix_Ln_distance;
   template <typename T> class matrix_X2_distance;
   template <typename T> class matrix_equal;
} /* namespace functors */
} /* namespace matrices */
} /* namespace math */
 
/*
 * Matrix implementation.
 */
namespace math {
namespace matrices {
/*
 * Imports.
 */
using concurrent::threads::synchronization::locks::auto_read_lock;
using concurrent::threads::synchronization::locks::auto_read_read_lock;
using concurrent::threads::synchronization::locks::auto_read_write_lock;
using concurrent::threads::synchronization::locks::auto_write_lock;
using concurrent::threads::synchronization::synchronizables::synchronizable;
using concurrent::threads::synchronization::synchronizables::unsynchronized;
using ::functors::comparable_functor;
using ::functors::compare_functors;
using io::serialization::serial_input_stream;
using io::serialization::serial_output_stream;
using io::serialization::serializer;
using io::serialization::serializers;
using io::streams::ios;
using io::streams::ostream;
using io::streams::ostringstream;
using lang::array;
using lang::exceptions::ex_index_out_of_bounds;
using lang::exceptions::ex_invalid_argument;
using math::complex;
using math::matrices::exceptions::ex_matrix_singular;
using math::matrices::exceptions::ex_matrix_dimension_mismatch;
using math::random::generators::rand_gen;
using math::random::generators::rand_gen_uniform;

/*
 * Import standard math functions.
 */
using math::cos;
using math::sin;
using math::tan;
using math::cosh;
using math::sinh;
using math::tanh;
using math::atan;
using math::asinh;
using math::exp;
using math::pow;
using math::atan2;
using math::abs;
using math::eps;

/*
 * Declare existence of matrix and complex matrix classes.
 */
template <typename T, typename Syn> class matrix;
template <typename T, typename Syn> class cmatrix;

/*
 * Declare prototypes for matrix template friend functions.
 */
template <typename T, typename Syn>
ostream& operator<<(ostream&, const matrix<T,Syn>&);

template <typename T, typename Syn>
matrix<T,Syn> operator+(const matrix<T,Syn>&, const T&);
template <typename T, typename Syn>
matrix<T,Syn> operator-(const matrix<T,Syn>&, const T&);
template <typename T, typename Syn>
matrix<T,Syn> operator*(const matrix<T,Syn>&, const T&);
template <typename T, typename Syn>
matrix<T,Syn> operator/(const matrix<T,Syn>&, const T&);

template <typename T, typename Syn>
matrix<T,Syn> operator+(const T&, const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<T,Syn> operator-(const T&, const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<T,Syn> operator*(const T&, const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<T,Syn> operator/(const T&, const matrix<T,Syn>&);

template <typename T, typename Syn>
matrix<T,Syn> operator+(const matrix<T,Syn>&, const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<T,Syn> operator-(const matrix<T,Syn>&, const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<T,Syn> operator*(const matrix<T,Syn>&, const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<T,Syn> operator/(const matrix<T,Syn>&, const matrix<T,Syn>&);

template <typename T, typename Syn>
matrix<bool,Syn> operator==(const matrix<T,Syn>&, const T&);
template <typename T, typename Syn>
matrix<bool,Syn> operator!=(const matrix<T,Syn>&, const T&);
template <typename T, typename Syn>
matrix<bool,Syn> operator< (const matrix<T,Syn>&, const T&);
template <typename T, typename Syn>
matrix<bool,Syn> operator> (const matrix<T,Syn>&, const T&);
template <typename T, typename Syn>
matrix<bool,Syn> operator<=(const matrix<T,Syn>&, const T&);
template <typename T, typename Syn>
matrix<bool,Syn> operator>=(const matrix<T,Syn>&, const T&);

template <typename T, typename Syn>
matrix<bool,Syn> operator==(const T&, const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator!=(const T&, const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator< (const T&, const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator> (const T&, const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator<=(const T&, const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator>=(const T&, const matrix<T,Syn>&);

template <typename T, typename Syn>
matrix<bool,Syn> operator==(const matrix<T,Syn>&, const complex<T>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator!=(const matrix<T,Syn>&, const complex<T>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator< (const matrix<T,Syn>&, const complex<T>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator> (const matrix<T,Syn>&, const complex<T>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator<=(const matrix<T,Syn>&, const complex<T>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator>=(const matrix<T,Syn>&, const complex<T>&);

template <typename T, typename Syn>
matrix<bool,Syn> operator==(const complex<T>&, const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator!=(const complex<T>&, const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator< (const complex<T>&, const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator> (const complex<T>&, const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator<=(const complex<T>&, const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator>=(const complex<T>&, const matrix<T,Syn>&);

template <typename T, typename Syn>
matrix<bool,Syn> operator==(const matrix<T,Syn>&, const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator!=(const matrix<T,Syn>&, const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator< (const matrix<T,Syn>&, const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator> (const matrix<T,Syn>&, const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator<=(const matrix<T,Syn>&, const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator>=(const matrix<T,Syn>&, const matrix<T,Syn>&);

template <typename T, typename Syn>
matrix<T,Syn> cos(const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<T,Syn> sin(const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<T,Syn> tan(const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<T,Syn> cosh(const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<T,Syn> sinh(const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<T,Syn> tanh(const matrix<T,Syn>&);

template <typename T, typename Syn>
cmatrix<T,Syn> acos(const matrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> asin(const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<T,Syn> atan(const matrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> acosh(const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<T,Syn> asinh(const matrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> atanh(const matrix<T,Syn>&);

template <typename T, typename Syn>
cmatrix<T,Syn> sqrt(const matrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> log(const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<T,Syn> exp(const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<T,Syn> pow(const matrix<T,Syn>&, int);
template <typename T, typename Syn>
cmatrix<T,Syn> pow(const matrix<T,Syn>&, const T&);
template <typename T, typename Syn>
cmatrix<T,Syn> pow(const matrix<T,Syn>&, const complex<T>&);

template <typename T, typename Syn>
matrix<T,Syn> atan2(const matrix<T,Syn>&, const matrix<T,Syn>&);

template <typename T, typename Syn>
matrix<T,Syn> abs(const matrix<T,Syn>&);

template <typename T, typename Syn>
matrix<T,Syn> conv(const matrix<T,Syn>&, const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<T,Syn> conv_crop(const matrix<T,Syn>&, const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<T,Syn> conv_crop_strict(const matrix<T,Syn>&, const matrix<T,Syn>&);

template <typename T, typename Syn>
T dot(const matrix<T,Syn>&, const matrix<T,Syn>&);

template <typename T, typename Syn>
matrix<T,Syn> prod(const matrix<T,Syn>&, const matrix<T,Syn>&);

template <typename T, typename Syn> T min(const matrix<T,Syn>&);
template <typename T, typename Syn> T max(const matrix<T,Syn>&);
template <typename T, typename Syn> T sum(const matrix<T,Syn>&);
template <typename T, typename Syn> T prod(const matrix<T,Syn>&);

template <typename T, typename Syn>
matrix<T,Syn> min(const matrix<T,Syn>&, unsigned long);
template <typename T, typename Syn>
matrix<T,Syn> max(const matrix<T,Syn>&, unsigned long);
template <typename T, typename Syn>
matrix<T,Syn> sum(const matrix<T,Syn>&, unsigned long);
template <typename T, typename Syn>
matrix<T,Syn> prod(const matrix<T,Syn>&, unsigned long);

template <typename T, typename Syn>
matrix<T,Syn> cumsum(const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<T,Syn> cumprod(const matrix<T,Syn>&);

template <typename T, typename Syn>
matrix<T,Syn> cumsum(const matrix<T,Syn>&, unsigned long);
template <typename T, typename Syn>
matrix<T,Syn> cumprod(const matrix<T,Syn>&, unsigned long);

template <typename T, typename Syn> T mean(const matrix<T,Syn>&);
template <typename T, typename Syn> T var(const matrix<T,Syn>&);

template <typename T, typename Syn>
matrix<T,Syn> mean(const matrix<T,Syn>&, unsigned long);
template <typename T, typename Syn>
matrix<T,Syn> var(const matrix<T,Syn>&, unsigned long);

template <typename T, typename Syn>
matrix<T,Syn> gradient(const matrix<T,Syn>&, unsigned long);
template <typename T, typename Syn>
matrix<T,Syn> gradient(const matrix<T,Syn>&, unsigned long, const T&);

template <typename T, typename Syn>
matrix<T,Syn> diag(const matrix<T,Syn>&);

template <typename T, typename Syn>
matrix<T,Syn> inv(const matrix<T,Syn>&, const T& = eps<T>());
template <typename T, typename Syn>
matrix<T,Syn> ref(const matrix<T,Syn>&, const T& = eps<T>());
template <typename T, typename Syn>
matrix<T,Syn> rref(const matrix<T,Syn>&, const T& = eps<T>());

template <typename T, typename Syn>
matrix<T,Syn> transpose(const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<T,Syn> shift_dimensions(
   const matrix<T,Syn>&, unsigned long);
template <typename T, typename Syn>
matrix<T,Syn> permute_dimensions(
   const matrix<T,Syn>&, const array<unsigned long>&);
template <typename T, typename Syn>
matrix<T,Syn> squeeze_dimensions(const matrix<T,Syn>&);

template <typename T, typename Syn>
matrix<T,Syn> resize(const matrix<T,Syn>&, unsigned long);
template <typename T, typename Syn>
matrix<T,Syn> resize(const matrix<T,Syn>&, unsigned long, const T&);

template <typename T, typename Syn>
matrix<T,Syn> resize(
   const matrix<T,Syn>&, unsigned long, unsigned long);
template <typename T, typename Syn>
matrix<T,Syn> resize(
   const matrix<T,Syn>&, unsigned long, unsigned long, const T&);

template <typename T, typename Syn>
matrix<T,Syn> resize(
   const matrix<T,Syn>&, const array<unsigned long>&);
template <typename T, typename Syn>
matrix<T,Syn> resize(
   const matrix<T,Syn>&, const array<unsigned long>&, const T&);

template <typename T, typename Syn>
matrix<T,Syn> vector(const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<T,Syn> reshape(const matrix<T,Syn>&, unsigned long, unsigned long);
template <typename T, typename Syn>
matrix<T,Syn> reshape(const matrix<T,Syn>&, const array<unsigned long>&);

template <typename T, typename Syn>
matrix<T,Syn> repmat(const matrix<T,Syn>&, unsigned long, unsigned long);
template <typename T, typename Syn>
matrix<T,Syn> repmat(const matrix<T,Syn>&, const array<unsigned long>&);

template <typename T, typename Syn>
matrix<T,Syn> vertcat(const matrix<T,Syn>&, const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<T,Syn> horzcat(const matrix<T,Syn>&, const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<T,Syn> concat(const matrix<T,Syn>&, const matrix<T,Syn>&, unsigned long);

/*
 * Base class containing code common to all matrix<T,Syn> templates.
 */
class matrix_base {
public:
   /*
    * Destructor.
    */
   virtual ~matrix_base() = 0;

protected:
   /*
    * Compute which dimension is shortest (return lower id in case of a tie).
    */
   static unsigned long dims_min(const array<unsigned long>&);

   /*
    * Compute which dimension is longest (return lower id in case of a tie).
    */
   static unsigned long dims_max(const array<unsigned long>&);

   /*
    * Compute size of matrix with given dimensions.
    */
   static unsigned long dims_matrix_size(const array<unsigned long>&);

   /*
    * Compute size of a step along the specified dimension (change in linear
    * index).
    */
   static unsigned long dims_size(const array<unsigned long>&, unsigned long);
   
   /*
    * Compute size of a step in each dimension (change in linear index).
    */
   static array<unsigned long> dims_sizes(const array<unsigned long>&);
   
   /*
    * Check if dimensions are equal.
    */
   static bool dims_equal(
      const array<unsigned long>&, 
      const array<unsigned long>&
   );

   /*
    * Assert that dimensions are equal.
    * Throw an exception (ex_matrix_dimension_mismatch) if they are not equal.
    */
   static void assert_dims_equal(
      const array<unsigned long>&, 
      const array<unsigned long>&
   );
   
   /*
    * Check if dimensions make matrices compatible for multiplcation.
    */
   static bool dims_mult_compatible(
      const array<unsigned long>&, 
      const array<unsigned long>&
   );

   /*
    * Compute dimensions resulting from matrix multiplication.
    * Throw an exception (ex_matrix_dimension_mismatch) if the dimensions 
    * are incompatible.
    */
   static array<unsigned long> dims_mult(
      const array<unsigned long>&, 
      const array<unsigned long>&
   );

   /*
    * Extend the dimensions array to at least the given number of dimensions
    * by appending dimensions of size one (as needed).
    */
   static array<unsigned long> extend_dims(
      const array<unsigned long>&,  /* dimensions */
      unsigned long                 /* minimum number of desired dimensions */
   );
 
   /*
    * Compute the range of indices in a block.
    */
   static array< array<unsigned long> > range_indices(
      const array<unsigned long>&,  /* start index along each dimension */
      const array<unsigned long>&   /* end index along each dimension */
   );

   /*
    * Compute the range of indices in a block.
    */
   static array< array<unsigned long> > range_indices(
      const array<unsigned long>&,  /* start index along each dimension */
      const array<unsigned long>&,  /* step along each dimension */
      const array<unsigned long>&   /* end index along each dimension */
   );

   /*
    * Compute linear index from two indices.
    */
   static unsigned long linear_index(
      const array<unsigned long>&,  /* dimensions */
      unsigned long,                /* indices */
      unsigned long
   );
   
   /*
    * Compute linear index from multiple indices.
    */
   static unsigned long linear_index(
      const array<unsigned long>&,  /* dimensions */
      const array<unsigned long>&   /* indices */
   );

   /*
    * Compute linear indices along the specified slice of the given dimenion.
    * These are the linear indices for all elements for which the coordinate 
    * along the given dimension is the given value (defaults to first slice).
    */
   static array<unsigned long> linear_indices_slice(
      const array<unsigned long>&,  /* dimensions */
      unsigned long,                /* dimension along which to slice */
      unsigned long = 0             /* coordinate of slice */
   );

   /*
    * Compute linear indices for a block of indices.
    */
   static array<unsigned long> linear_indices(
      const array<unsigned long>&,  /* dimensions */
      const array<unsigned long>&,  /* start index along each dimension */
      const array<unsigned long>&   /* end index along each dimension */
   );

   /*
    * Compute linear indices for a block of indices.
    */
   static array<unsigned long> linear_indices(
      const array<unsigned long>&,  /* dimensions */
      const array<unsigned long>&,  /* start index along each dimension */
      const array<unsigned long>&,  /* step along each dimension */
      const array<unsigned long>&   /* end index along each dimension */
   );

   /*
    * Compute linear indices for an arbitrary range of indices.
    */
   static array<unsigned long> linear_indices(
      const array<unsigned long>&,           /* dimensions */
      const array< array<unsigned long> >&   /* indices along each dimension */
   );
   
   /*
    * Compute linear indices for an arbitrary range of indices and an
    * arbitrary matrix memory layout.  Use the specified step sizes along 
    * each dimension rather than assuming the standard memory layout.
    */
   static array<unsigned long> linear_indices_custom(
      const array<unsigned long>&,           /* step sizes */
      const array< array<unsigned long> >&   /* indices along each dimension */
   );

   /*
    * Computing relative linear offsets of neighbors from the center of a
    * block of 3^d elements for a d-dimensional matrix of the given dimensions.
    *
    * Return an array of 3^d offsets (including the center of the block).
    * Note that some offsets will be negative.
    */
   static array<long> linear_neighbor_offsets(
      const array<unsigned long>&   /* dimensions */
   );
};

/*
 * Matrix class.
 */
template <typename T = double, typename Syn = unsynchronized>
class matrix : public array<T,Syn>, 
               protected matrix_base {
public:
   /*
    * Friend classes.
    */
   template <typename U, typename S> friend class matrix;
   template <typename U, typename S> friend class cmatrix;

   /*
    * Friend libraries.
    */
   friend class math::libraries::lib_image;
   friend class math::libraries::lib_matrix;
   friend class math::libraries::lib_signal;

   /*
    * Friend functors.
    */
   friend class math::matrices::functors::matrix_L1_distance<T>;
   friend class math::matrices::functors::matrix_L2_distance<T>;
   friend class math::matrices::functors::matrix_Ln_distance<T>;
   friend class math::matrices::functors::matrix_X2_distance<T>;
   friend class math::matrices::functors::matrix_equal<T>;

   /*
    * Constructor.
    * Return the empty (0 x 0) matrix.
    */
   matrix();
   
   /*
    * Constructor for single element (1 x 1) matrix.
    */
   explicit matrix(const T&);
   
   /*
    * Constructors for two-dimensional M x N matrices.
    */
   explicit matrix(
      unsigned long,                /* M */
      unsigned long                 /* N */
   );
   
   explicit matrix(
      unsigned long,                /* M */
      unsigned long,                /* N */
      const T&                      /* element initialization value */
   );
   
   /*
    * Constructors for multi-dimensional matrices.
    */
   explicit matrix(
      const array<unsigned long>&   /* dimensions */
   );
      
   explicit matrix(
      const array<unsigned long>&,  /* dimensions */
      const T&                      /* element initialization value */
   );

   /*
    * Copy constructors.
    */
   matrix(const matrix<T,Syn>&);
   template <typename U, typename S> explicit matrix(const matrix<U,S>&);
   
   /*
    * Destructor.
    */
   virtual ~matrix();
   
   /*
    * Empty matrix (0 x 0).
    */
   static matrix<T,Syn> empty();
   
   /*
    * Identity matrix.
    * Create and return an N x N identity matrix.
    */
   static matrix<T,Syn> eye(unsigned long /* N */);
   
   /*
    * Zeros matrix.
    * Create and return an M x N matrix of zeros.
    */
   static matrix<T,Syn> zeros(
      unsigned long,                /* M */
      unsigned long                 /* N */
   );

   /*
    * Zeros matrix.
    * Create and return a matrix of zeros with the given dimensions.
    */
   static matrix<T,Syn> zeros(
      const array<unsigned long>&   /* dimensions */
   );
   
   /*
    * Ones matrix.
    * Create and return an M x N matrix of ones.
    */
   static matrix<T,Syn> ones(
      unsigned long,                /* M */
      unsigned long                 /* N */
   );
   
   /*
    * Ones matrix.
    * Create and return a matrix of ones with the given dimensions.
    */
   static matrix<T,Syn> ones(
      const array<unsigned long>&   /* dimensions */
   );

   /*
    * Random matrix.
    * Create and return an M x N matrix of random numbers drawn from the
    * given random number generator.  If no generator is specified, then
    * the matrix elements are uniformly distributed random numbers in the 
    * interval [0, 1] (inclusive).
    */ 
   static matrix<T,Syn> random(
      unsigned long,                /* M */
      unsigned long                 /* N */
   );

   static matrix<T,Syn> random(
      unsigned long,                /* M */
      unsigned long,                /* N */
      rand_gen<T>&                  /* generator */
   );
  
   /*
    * Random matrix.
    * Create and return a matrix with the given dimensions in which each 
    * element is a random number drawn from the given random number generator.
    * If no generator is specified, then the matrix elements are uniformly 
    * distributed random numbers in the interval [0, 1] (inclusive).
    */
   static matrix<T,Syn> random(
      const array<unsigned long>&   /* dimensions */
   );

   static matrix<T,Syn> random(
      const array<unsigned long>&,  /* dimensions */
      rand_gen<T>&                  /* generator */
   );

   /*
    * Ramp.
    * Create and return a vector containing values from start to end at an 
    * increment of step.  The default step is one.
    */
   static matrix<T,Syn> ramp(
      const T&,                     /* start */
      const T&                      /* end */
   );
   
   static matrix<T,Syn> ramp(
      const T&,                     /* start */
      const T&,                     /* step */
      const T&                      /* end */
   );
   
   /*
    * Diagonal matrix.
    * Create and return a d-dimensional square matrix given its diagonal.
    */
   static matrix<T,Syn> diagonal(
      const matrix<T,Syn>&,         /* diagonal entries */
      unsigned long = 2             /* d (dimensionality) */
   );

   /*
    * Serialize.
    */
   void serialize(
      serial_output_stream&,
      const serializer<T>& = serializers<T>::s_default()
   ) const;

   /*
    * Deserialize.
    */
   static auto_ptr< matrix<T,Syn> > deserialize(
      serial_input_stream&,
      const serializer<T>& = serializers<T>::s_default()
   );

   /*
    * Element reference (using two indices).
    */
   T& operator()(unsigned long, unsigned long);
   const T& operator()(unsigned long, unsigned long) const;
    
   /*
    * Element reference (using multiple indices).
    */
   T& operator()(const array<unsigned long>&);
   const T& operator()(const array<unsigned long>&) const;
   
   /*
    * Direct element access.
    * Return a pointer to the memory containing the matrix elements.
    * This should be used with extreme caution.
    */
   T* data();
   const T* data() const;

   /*
    * Subarray.
    * Return vector of elements in the specified linear index range.
    */
   matrix<T,Syn> subarray(
      unsigned long,                /* start index */ 
      unsigned long                 /* end index   */
   ) const;
   
   matrix<T,Syn> subarray(
      unsigned long,                /* start index */ 
      unsigned long,                /* step size   */
      unsigned long                 /* end index   */
   ) const;
   
   matrix<T,Syn> subarray(
      const array<unsigned long>&   /* indices of desired elements */
   ) const;
  
   /*
    * Submatrix.
    * Return matrix of elements in the specified range.
    */
   matrix<T,Syn> submatrix(
      const array<unsigned long>&,  /* start index along each dimension */
      const array<unsigned long>&   /* end index along each dimension */
   ) const;
   
   matrix<T,Syn> submatrix(
      const array<unsigned long>&,  /* start index along each dimension */
      const array<unsigned long>&,  /* step along each dimension */
      const array<unsigned long>&   /* end index along each dimension */
   ) const;
   
   matrix<T,Syn> submatrix(
      const array< array<unsigned long> >&   /* indices along each dimension */
   ) const;

   /*
    * Subarray assignment.
    * Assign value(s) to elements in the specified linear index range.
    */
   matrix<T,Syn>& subassign(
      unsigned long,                /* start index */ 
      unsigned long,                /* end index   */
      const T&                      /* value */
   );

   matrix<T,Syn>& subassign(
      unsigned long,                /* start index */ 
      unsigned long,                /* end index   */
      const matrix<T,Syn>&          /* values */
   );

   matrix<T,Syn>& subassign(
      unsigned long,                /* start index */ 
      unsigned long,                /* step size   */
      unsigned long,                /* end index   */
      const T&                      /* value */
   );

   matrix<T,Syn>& subassign(
      unsigned long,                /* start index */ 
      unsigned long,                /* step size   */
      unsigned long,                /* end index   */
      const matrix<T,Syn>&          /* values */
   );

   matrix<T,Syn>& subassign(
      const array<unsigned long>&,  /* indices of desired elements */
      const T&                      /* value */
   );

   matrix<T,Syn>& subassign(
      const array<unsigned long>&,  /* indices of desired elements */
      const matrix<T,Syn>&          /* values */
   );

   /*
    * Submatrix assignment.
    * Assign value(s) to elements in the specified range.
    */
   matrix<T,Syn>& subassign(
      const array<unsigned long>&,  /* start index along each dimension */
      const array<unsigned long>&,  /* end index along each dimension */
      const T&                      /* value */
   );

   matrix<T,Syn>& subassign(
      const array<unsigned long>&,  /* start index along each dimension */
      const array<unsigned long>&,  /* end index along each dimension */
      const matrix<T,Syn>&          /* values */
   );
 
   matrix<T,Syn>& subassign(
      const array<unsigned long>&,  /* start index along each dimension */
      const array<unsigned long>&,  /* step along each dimension */
      const array<unsigned long>&,  /* end index along each dimension */
      const T&                      /* value */
   );
     
   matrix<T,Syn>& subassign(
      const array<unsigned long>&,  /* start index along each dimension */
      const array<unsigned long>&,  /* step along each dimension */
      const array<unsigned long>&,  /* end index along each dimension */
      const matrix<T,Syn>&          /* values */
   );
    
   matrix<T,Syn>& subassign(
      const array< array<unsigned long> >&,  /* indices along each dimension */
      const T&                               /* value */
   );
  
   matrix<T,Syn>& subassign(
      const array< array<unsigned long> >&,  /* indices along each dimension */
      const matrix<T,Syn>&                   /* values */
   );
   
   /*
    * Fill with value.
    */
   matrix<T,Syn>& fill(const T&);

   /*
    * Formatted output to stream.
    */
   friend ostream& operator<< <T,Syn>(ostream&, const matrix<T,Syn>&);
   
   /*
    * Assignment operators: matrix-real.
    */
   matrix<T,Syn>& operator=(const T&);
   matrix<T,Syn>& operator+=(const T&);
   matrix<T,Syn>& operator-=(const T&);
   matrix<T,Syn>& operator*=(const T&);
   matrix<T,Syn>& operator/=(const T&);
   
   /*
    * Assignment operators: matrix-matrix.
    */
   matrix<T,Syn>& operator=(const matrix<T,Syn>&);
   
   template <typename U, typename S> matrix<T,Syn>&
      operator=(const matrix<U,S>&);
   template <typename U, typename S> matrix<T,Syn>&
      operator+=(const matrix<U,S>&); 
   template <typename U, typename S> matrix<T,Syn>&
      operator-=(const matrix<U,S>&);
   template <typename U, typename S> matrix<T,Syn>&
      operator*=(const matrix<U,S>&);
   template <typename U, typename S> matrix<T,Syn>&
      operator/=(const matrix<U,S>&);
   
   /*
    * Binary matrix-real operators.
    */
   friend matrix<T,Syn> operator+ <T,Syn>(const matrix<T,Syn>&, const T&);
   friend matrix<T,Syn> operator- <T,Syn>(const matrix<T,Syn>&, const T&);
   friend matrix<T,Syn> operator* <T,Syn>(const matrix<T,Syn>&, const T&);
   friend matrix<T,Syn> operator/ <T,Syn>(const matrix<T,Syn>&, const T&);

   friend matrix<T,Syn> operator+ <T,Syn>(const T&, const matrix<T,Syn>&);
   friend matrix<T,Syn> operator- <T,Syn>(const T&, const matrix<T,Syn>&);
   friend matrix<T,Syn> operator* <T,Syn>(const T&, const matrix<T,Syn>&);
   friend matrix<T,Syn> operator/ <T,Syn>(const T&, const matrix<T,Syn>&);
   
   /*
    * Binary matrix-matrix operators.
    */
   friend matrix<T,Syn>
      operator+ <T,Syn>(const matrix<T,Syn>&, const matrix<T,Syn>&);
   friend matrix<T,Syn>
      operator- <T,Syn>(const matrix<T,Syn>&, const matrix<T,Syn>&);
   friend matrix<T,Syn>
      operator* <T,Syn>(const matrix<T,Syn>&, const matrix<T,Syn>&);
   friend matrix<T,Syn>
      operator/ <T,Syn>(const matrix<T,Syn>&, const matrix<T,Syn>&);

   /*
    * Binary matrix-real comparators.
    */
   friend matrix<bool,Syn> operator== <T,Syn>(const matrix<T,Syn>&, const T&);
   friend matrix<bool,Syn> operator!= <T,Syn>(const matrix<T,Syn>&, const T&);
   friend matrix<bool,Syn> operator<  <T,Syn>(const matrix<T,Syn>&, const T&);
   friend matrix<bool,Syn> operator>  <T,Syn>(const matrix<T,Syn>&, const T&);
   friend matrix<bool,Syn> operator<= <T,Syn>(const matrix<T,Syn>&, const T&);
   friend matrix<bool,Syn> operator>= <T,Syn>(const matrix<T,Syn>&, const T&);
   
   friend matrix<bool,Syn> operator== <T,Syn>(const T&, const matrix<T,Syn>&);
   friend matrix<bool,Syn> operator!= <T,Syn>(const T&, const matrix<T,Syn>&);
   friend matrix<bool,Syn> operator<  <T,Syn>(const T&, const matrix<T,Syn>&);
   friend matrix<bool,Syn> operator>  <T,Syn>(const T&, const matrix<T,Syn>&);
   friend matrix<bool,Syn> operator<= <T,Syn>(const T&, const matrix<T,Syn>&);
   friend matrix<bool,Syn> operator>= <T,Syn>(const T&, const matrix<T,Syn>&);
      
   /*
    * Binary matrix-complex comparators.
    */
   friend matrix<bool,Syn>
      operator== <T,Syn>(const matrix<T,Syn>&, const complex<T>&);
   friend matrix<bool,Syn>
      operator!= <T,Syn>(const matrix<T,Syn>&, const complex<T>&);
   friend matrix<bool,Syn>
      operator<  <T,Syn>(const matrix<T,Syn>&, const complex<T>&);
   friend matrix<bool,Syn>
      operator>  <T,Syn>(const matrix<T,Syn>&, const complex<T>&);
   friend matrix<bool,Syn>
      operator<= <T,Syn>(const matrix<T,Syn>&, const complex<T>&);
   friend matrix<bool,Syn>
      operator>= <T,Syn>(const matrix<T,Syn>&, const complex<T>&);
   
   friend matrix<bool,Syn>
      operator== <T,Syn>(const complex<T>&, const matrix<T,Syn>&);
   friend matrix<bool,Syn>
      operator!= <T,Syn>(const complex<T>&, const matrix<T,Syn>&);
   friend matrix<bool,Syn>
      operator<  <T,Syn>(const complex<T>&, const matrix<T,Syn>&);
   friend matrix<bool,Syn>
      operator>  <T,Syn>(const complex<T>&, const matrix<T,Syn>&);
   friend matrix<bool,Syn>
      operator<= <T,Syn>(const complex<T>&, const matrix<T,Syn>&);
   friend matrix<bool,Syn>
      operator>= <T,Syn>(const complex<T>&, const matrix<T,Syn>&);
      
   /*
    * Binary matrix-matrix comparators.
    */
   friend matrix<bool,Syn>
      operator== <T,Syn>(const matrix<T,Syn>&, const matrix<T,Syn>&);
   friend matrix<bool,Syn>
      operator!= <T,Syn>(const matrix<T,Syn>&, const matrix<T,Syn>&);
   friend matrix<bool,Syn>
      operator<  <T,Syn>(const matrix<T,Syn>&, const matrix<T,Syn>&);
   friend matrix<bool,Syn>
      operator>  <T,Syn>(const matrix<T,Syn>&, const matrix<T,Syn>&);
   friend matrix<bool,Syn>
      operator<= <T,Syn>(const matrix<T,Syn>&, const matrix<T,Syn>&);
   friend matrix<bool,Syn>
      operator>= <T,Syn>(const matrix<T,Syn>&, const matrix<T,Syn>&);
 
   /*
    * Unary arithmetic operators.
    */
   matrix<T,Syn> operator+() const;
   matrix<T,Syn> operator-() const;
   matrix<T,Syn>& operator++();
   matrix<T,Syn>& operator--();
   matrix<T,Syn> operator++(int);
   matrix<T,Syn> operator--(int);
   
   /*
    * Unary element logical inversion operator.
    */
   matrix<bool,Syn> operator!() const;
   
   /*
    * Unary element multiplicative inversion operator.
    */
   matrix<T,Syn> operator~() const;
  
   /*
    * Transcendentals on elements.
    */
   friend  matrix<T,Syn> cos<T,Syn>(const matrix<T,Syn>&);
   friend  matrix<T,Syn> sin<T,Syn>(const matrix<T,Syn>&);
   friend  matrix<T,Syn> tan<T,Syn>(const matrix<T,Syn>&);
   friend  matrix<T,Syn> cosh<T,Syn>(const matrix<T,Syn>&);
   friend  matrix<T,Syn> sinh<T,Syn>(const matrix<T,Syn>&);
   friend  matrix<T,Syn> tanh<T,Syn>(const matrix<T,Syn>&);
   
   friend cmatrix<T,Syn> acos<T,Syn>(const matrix<T,Syn>&);
   friend cmatrix<T,Syn> asin<T,Syn>(const matrix<T,Syn>&);
   friend  matrix<T,Syn> atan<T,Syn>(const matrix<T,Syn>&);
   friend cmatrix<T,Syn> acosh<T,Syn>(const matrix<T,Syn>&);
   friend  matrix<T,Syn> asinh<T,Syn>(const matrix<T,Syn>&);
   friend cmatrix<T,Syn> atanh<T,Syn>(const matrix<T,Syn>&);

   friend cmatrix<T,Syn> sqrt<T,Syn>(const matrix<T,Syn>&);
   friend cmatrix<T,Syn> log<T,Syn>(const matrix<T,Syn>&);
   friend  matrix<T,Syn> exp<T,Syn>(const matrix<T,Syn>&);
   friend  matrix<T,Syn> pow<T,Syn>(const matrix<T,Syn>&, int);
   friend cmatrix<T,Syn> pow<T,Syn>(const matrix<T,Syn>&, const T&);
   friend cmatrix<T,Syn> pow<T,Syn>(const matrix<T,Syn>&, const complex<T>&);
   /*
    * Two-parameter arctangent of the corresponding elements.
    */
   friend matrix<T,Syn> atan2<T,Syn>(
      const matrix<T,Syn>& /* y */, const matrix<T,Syn>& /* x */
   );
   
   /*
    * Absolute value of elements.
    */
   friend matrix<T,Syn> abs<T,Syn>(const matrix<T,Syn>&);
   
   /*
    * Convolution.
    */
   friend matrix<T,Syn>
      conv<T,Syn>(const matrix<T,Syn>&, const matrix<T,Syn>&);

   /*
    * Convolution.
    * Crop the result to be no larger than the left input matrix.
    */
   friend matrix<T,Syn>
      conv_crop<T,Syn>(const matrix<T,Syn>&, const matrix<T,Syn>&);

   /*
    * Convolution.
    * Crop the result to be no larger than the left input matrix and to
    * include only the central portion for which no zero padding was required.
    */
   friend matrix<T,Syn>
      conv_crop_strict<T,Syn>(const matrix<T,Syn>&, const matrix<T,Syn>&);
   
   /*
    * Dot product.
    */
   friend T dot<T,Syn>(const matrix<T,Syn>&, const matrix<T,Syn>&);
  
   /*
    * Element product.
    * Compute the product of corresponding matrix elements.
    */
   friend matrix<T,Syn> prod<T,Syn>(const matrix<T,Syn>&, const matrix<T,Syn>&);
   
   /*
    * Minimum, maximum, sum, product of all matrix elements.
    */
   friend T min<T,Syn>(const matrix<T,Syn>&);
   friend T max<T,Syn>(const matrix<T,Syn>&);
   friend T sum<T,Syn>(const matrix<T,Syn>&);
   friend T prod<T,Syn>(const matrix<T,Syn>&);
    
   /*
    * Minimum, maximum, sum, product along specified dimension.
    */
   friend matrix<T,Syn> min<T,Syn>(const matrix<T,Syn>&, unsigned long);
   friend matrix<T,Syn> max<T,Syn>(const matrix<T,Syn>&, unsigned long);
   friend matrix<T,Syn> sum<T,Syn>(const matrix<T,Syn>&, unsigned long);
   friend matrix<T,Syn> prod<T,Syn>(const matrix<T,Syn>&, unsigned long);
 
   /*
    * Cumulative sum, product of all matrix elements.
    */
   friend matrix<T,Syn> cumsum<T,Syn>(const matrix<T,Syn>&);
   friend matrix<T,Syn> cumprod<T,Syn>(const matrix<T,Syn>&);

   /*
    * Cumulative sum, product along specified dimension.
    */
   friend matrix<T,Syn> cumsum<T,Syn>(const matrix<T,Syn>&, unsigned long);
   friend matrix<T,Syn> cumprod<T,Syn>(const matrix<T,Syn>&, unsigned long);

   /*
    * Mean, variance of all matrix elements.
    */
   friend T mean<T,Syn>(const matrix<T,Syn>&);
   friend T var<T,Syn>(const matrix<T,Syn>&);

   /*
    * Mean, variance along specified dimension.
    */
   friend matrix<T,Syn> mean<T,Syn>(const matrix<T,Syn>&, unsigned long);
   friend matrix<T,Syn> var<T,Syn>(const matrix<T,Syn>&, unsigned long);
   
   /*
    * Gradient along specified dimension.
    * Optionally specify the spacing (defaults to unit spacing) of elements
    * along the dimension.
    */
   friend matrix<T,Syn>
      gradient<T,Syn>(const matrix<T,Syn>&, unsigned long);
   friend matrix<T,Syn>
      gradient<T,Syn>(const matrix<T,Syn>&, unsigned long, const T&);
   
   /*
    * Logical AND of elements (true iff all element(s) are nonzero).
    */
   bool all() const;
   
   /*
    * Logical OR of elements (true iff at least one element is nonzero).
    */
   bool some() const;
   
   /*
    * Logical NOR of elements (true iff all element(s) are zero).
    */
   bool none() const;
         
   /*
    * Exactly one element true (true iff exactly one element is nonzero).
    */
   bool one() const;

   /*
    * Find locations of nonzero elements.
    */
   array<unsigned long> find() const;
   
   /*
    * Find locations of elements with the specified value.
    */
   array<unsigned long> find_value(const T&) const;

   /*
    * Diagonal entries (returns diagonal of matrix in a vector).
    */
   friend matrix<T,Syn> diag<T,Syn>(const matrix<T,Syn>&);
   matrix<T,Syn> diag() const;

   /*
    * Inverse.
    * Optionally specify the numeric tolerance.
    */
   friend matrix<T,Syn> inv<T,Syn>(const matrix<T,Syn>&, const T& /* tol */);
   matrix<T,Syn>& inv(const T& = eps<T>() /* tol */);
   
   /*
    * Row echelon form.
    * Optionally specify the numeric tolerance.
    */
   friend matrix<T,Syn> ref<T,Syn>(const matrix<T,Syn>&, const T& /* tol */);
   matrix<T,Syn>& ref(const T& = eps<T>() /* tol */);
   
   /*
    * Reduced row echelon form.
    * Optionally specify the numeric tolerance.
    */
   friend matrix<T,Syn> rref<T,Syn>(const matrix<T,Syn>&, const T& /* tol */);
   matrix<T,Syn>& rref(const T& = eps<T>() /* tol */);
   
   /*
    * Determinant.
    * Optionally specify the numeric tolerance.
    */
   T det(const T& = eps<T>() /* tol */) const;
   
   /*
    * Rank.
    * Optionally specify the numeric tolerance.
    */
   unsigned long rank(const T& = eps<T>() /* tol */) const;

   /*
    * Get dimensionality of matrix.
    */
   unsigned long dimensionality() const;

   /*
    * Get dimensions of matrix.
    */
   array<unsigned long> dimensions() const;
   
   /*
    * Transpose.
    * Reverse dimension order.
    * The input matrix is regarded as being at least two-dimensional.
    */
   friend matrix<T,Syn> transpose<T,Syn>(const matrix<T,Syn>&);
   matrix<T,Syn>& transpose();
  
   /*
    * Dimension shifting.
    * Shift dimensions of the matrix to the left, wrapping the leading
    * dimensions around to the right.
    */
   friend matrix<T,Syn> shift_dimensions<T,Syn>(
      const matrix<T,Syn>&,
      unsigned long                 /* number of dimensions to shift */
   );

   matrix<T,Syn>& shift_dimensions(
      unsigned long                 /* number of dimensions to shift */
   );
   
   /*
    * Dimension permutation.
    */
   friend matrix<T,Syn> permute_dimensions<T,Syn>(
      const matrix<T,Syn>&,
      const array<unsigned long>&   /* ordering of dimensions */
   );

   matrix<T,Syn>& permute_dimensions(
      const array<unsigned long>&   /* ordering of dimensions */
   );

   /*
    * Singleton dimension elimination.
    * Remove singleton dimensions, leaving the matrix elements unchanged.
    */
   friend matrix<T,Syn> squeeze_dimensions<T,Syn>(const matrix<T,Syn>&);
   matrix<T,Syn>& squeeze_dimensions();
    
   /*
    * Get total size (number of matrix elements).
    */
   using array<T,Syn>::size;
   
   /*
    * Get size along specific dimension.
    */
   unsigned long size(unsigned long) const;
   
   /*
    * Resize (vector).
    * Return the resized vector.
    */
   friend matrix<T,Syn> resize<T,Syn>(
      const matrix<T,Syn>&,
      unsigned long                 /* length */
   );

   friend matrix<T,Syn> resize<T,Syn>(
      const matrix<T,Syn>&,
      unsigned long,                /* length */
      const T&                      /* element init value (if enlarging) */
   );
   
   /*
    * Resize (2D matrix).
    * Return the resized M x N matrix.
    */
   friend matrix<T,Syn> resize<T,Syn>(
      const matrix<T,Syn>&,
      unsigned long,                /* M */
      unsigned long                 /* N */
   );
   
   friend matrix<T,Syn> resize<T,Syn>(
      const matrix<T,Syn>&,
      unsigned long,                /* M */
      unsigned long,                /* N */
      const T&                      /* element init value (if enlarging) */
   );
   
   /*
    * Resize (multi-dimensional matrix).
    * Return the resized matrix.
    * The dimensionality of the matrix must not change.
    */
   friend matrix<T,Syn> resize<T,Syn>(
      const matrix<T,Syn>&,
      const array<unsigned long>&   /* dimensions */
   );
   
   friend matrix<T,Syn> resize<T,Syn>(
      const matrix<T,Syn>&,
      const array<unsigned long>&,  /* dimensions */
      const T&                      /* element init value (if enlarging) */
   );

   /*
    * Resize (vector).
    * The original vector is resized to the specified length.
    */
   matrix<T,Syn>& resize(
      unsigned long                 /* length */
   );

   matrix<T,Syn>& resize(
      unsigned long,                /* length */
      const T&                      /* element init value (if enlarging) */
   );

   /*
    * Resize (2D matrix).
    * The original 2D matrix is resized to be M x N.
    */
   matrix<T,Syn>& resize(
      unsigned long,                /* M */
      unsigned long                 /* N */
   );
   
   matrix<T,Syn>& resize(
      unsigned long,                /* M */
      unsigned long,                /* N */
      const T&                      /* element init value (if enlarging) */
   );
   
   /*
    * Resize (multi-dimensional matrix).
    * The original matrix is resized to the given dimensions.
    * The dimensionality of the matrix must not change.
    */
   matrix<T,Syn>& resize(
      const array<unsigned long>&   /* dimensions */
   );
   
   matrix<T,Syn>& resize(
      const array<unsigned long>&,  /* dimensions */
      const T&                      /* element init value (if enlarging) */
   );
     
   /*
    * Reshape into a vector.
    */
   friend matrix<T,Syn> vector<T,Syn>(const matrix<T,Syn>&);
   
   /*
    * Reshape into 2D M x N matrix.
    */
   friend matrix<T,Syn> reshape<T,Syn>(
      const matrix<T,Syn>&,
      unsigned long,                /* M */
      unsigned long                 /* N */
   );
   
   /*
    * Reshape into multi-dimensional matrix.
    */
   friend matrix<T,Syn> reshape<T,Syn>(
      const matrix<T,Syn>&,
      const array<unsigned long>&   /* dimensions */
   );
   
   /*
    * Reshape (in place) into a vector.
    */
   matrix<T,Syn>& vector();
  
   /*
    * Reshape (in place) into 2D M x N matrix.
    */
   matrix<T,Syn>& reshape(
      unsigned long,                /* M */
      unsigned long                 /* N */
   );
   
   /*
    * Reshape (in place) into multi-dimensional matrix.
    */
   matrix<T,Syn>& reshape(
      const array<unsigned long>&   /* dimensions */
   );
   
   /*
    * Replicate matrix.
    */
   friend matrix<T,Syn> repmat<T,Syn>(
      const matrix<T,Syn>&,
      unsigned long,                /* M */
      unsigned long                 /* N */
   );
      
   friend matrix<T,Syn> repmat<T,Syn>(
      const matrix<T,Syn>&, 
      const array<unsigned long>&   /* dimensions */
   );
   
   /*
    * Concatenate matrices.
    * Throw an exception (ex_matrix_dimension_mismatch) if the matrices
    * cannot be concatenated along the specified dimension.
    */
   friend matrix<T,Syn> vertcat<T,Syn>(
      const matrix<T,Syn>&, const matrix<T,Syn>&
   );
   friend matrix<T,Syn> horzcat<T,Syn>(
      const matrix<T,Syn>&, const matrix<T,Syn>&
   );
   friend matrix<T,Syn> concat<T,Syn>(
      const matrix<T,Syn>&, 
      const matrix<T,Syn>&, 
      unsigned long           /* dimension along which to concatenate */
   );
   
   /*
    * Reverse element order (treating the matrix as a linear array).
    */
   matrix<T,Syn>& reverse();

   /*
    * Reverse element order along the specified dimension.
    */
   matrix<T,Syn>& reverse(unsigned long /* dimension */);
    
   /*
    * Randomly permute the elements of the matrix, treating it as a linear
    * array.  Return an index array mapping resulting position --> original
    * position.
    */
   using array<T,Syn>::randperm;

   /*
    * Randomly permute along the specified dimension.
    * Return an index matrix mapping resulting position --> original position.
    */
   matrix<unsigned long> randperm(
      unsigned long  /* dimension */
   );
 
   matrix<unsigned long> randperm(
      unsigned long, /* dimension */
      rand_gen<>&    /* generator */
   );

   /*
    * Sort the matrix, treating it as a linear array.
    * Elements are sorted in ascending order according to the given comparison
    * functor.
    *
    * Sorting uses an O(n + (n/p)*log(n/p)) expected time algorithm, where
    * n is the array size and p is the number of available processors.
    */
   using array<T,Syn>::sort;

   /*
    * Sort the matrix, treating it as a linear array.
    * Return an index array mapping sorted position --> original position.
    */
   using array<T,Syn>::sort_idx;

   /*
    * Sort along specified dimension.
    * Elements are sorted in ascending order according to the given comparison
    * functor.
    */
   void sort(
      unsigned long, /* dimension */
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   );
   
   /*
    * Sort along specified dimension.
    * Return an index matrix mapping sorted position --> original position.
    */
   matrix<unsigned long> sort_idx(
      unsigned long, /* dimension */
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   );

   /*
    * Remove duplicate elements (as defined by the given comparison functor)
    * and arrange the remaining unique elements in sorted order as a vector.
    */
   void unique(
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   );

   /*
    * Remove duplicate elements (as defined by the given comparison functor)
    * and arrange the remaining unique elements in sorted order as a vector.
    * In addition, return an index array containing the original positions
    * of the unique elements.
    */
   array<unsigned long> unique_idx(
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   );

protected:
   /************************************************************************
    * Matrix data.
    ************************************************************************/

   array<unsigned long> _dims;   /* dimensions of matrix */

   /************************************************************************
    * Matrix helper functions/methods.
    * Note: These functions do not lock their matrix arguments.
    ************************************************************************/

   /*
    * Copy the contents of a matrix.
    */
   static matrix<T,Syn> copy(const matrix<T,Syn>&);

   /*
    * Copy the portion of the source matrix that overlaps the destination
    * matrix into the destination matrix.
    */
   static void copy_overlap(
      const matrix<T,Syn>&,   /* source */
      matrix<T,Syn>&          /* destination */
   );
      
   /*
    * Swap the contents of two matrices.
    */
   static void swap(matrix<T,Syn>&, matrix<T,Syn>&);
   
   /*
    * Convert an array<T,Syn> object into a matrix<T,Syn> object.
    * The original array<T,Syn> object is destroyed in the process.
    */
   static matrix<T,Syn> to_matrix(array<T,Syn>&);

   /*
    * Output classname to stream.
    */
   virtual ostream& output_classname(ostream&) const;

   /*
    * Output matrix element at given index to stream.
    */
   virtual ostream& output_element(ostream&, unsigned long) const;
   
   /*
    * Computation helper functions.
    */
   static matrix<T,Syn> compute_conv_1D(
      const matrix<T,Syn>&, const matrix<T,Syn>&, bool, bool);
   static matrix<T,Syn> compute_conv_2D(
      const matrix<T,Syn>&, const matrix<T,Syn>&,
      const array<unsigned long>&, const array<unsigned long>&, bool, bool);
   static matrix<T,Syn> compute_conv(
      const matrix<T,Syn>&, const matrix<T,Syn>&, bool, bool);
      
   static matrix<T,Syn> compute_inv(const matrix<T,Syn>&, const T&);
   
   static matrix<T,Syn> compute_transpose(const matrix<T,Syn>&);
   static matrix<T,Syn> compute_shift_dimensions(
      const matrix<T,Syn>&, unsigned long);
   static matrix<T,Syn> compute_permute_dimensions(
      const matrix<T,Syn>&, const array<unsigned long>&);
      
   static matrix<T,Syn> compute_resize(
      const matrix<T,Syn>&, unsigned long);
   static matrix<T,Syn> compute_resize(
      const matrix<T,Syn>&, unsigned long, const T&);

   static matrix<T,Syn> compute_resize(
      const matrix<T,Syn>&, unsigned long, unsigned long);
   static matrix<T,Syn> compute_resize(
      const matrix<T,Syn>&, unsigned long, unsigned long, const T&);

   static matrix<T,Syn> compute_resize(
      const matrix<T,Syn>&, const array<unsigned long>&);
   static matrix<T,Syn> compute_resize(
      const matrix<T,Syn>&, const array<unsigned long>&, const T&);

   /*
    * Row reduce a 2D matrix (in place).
    *
    * Put the matrix in either row echelon form or reduced row echelon form.
    * Return the factor by which the determinant has been scaled.
    */
   T row_reduce(
      const T&,      /* numeric tolerance */
      bool = true    /* reduced row echelon form? */ 
   );
};

/***************************************************************************
 * Constructors and destructor.
 ***************************************************************************/

/*
 * Constructor.
 * Return the empty (0 x 0) matrix.
 */
template <typename T, typename Syn>
matrix<T,Syn>::matrix()
 : array<T,Syn>(),
   _dims(2)
{ }

/*
 * Constructor for single element (1 x 1) matrix.
 */
template <typename T, typename Syn>
matrix<T,Syn>::matrix(const T& t)
 : array<T,Syn>(1, t),
   _dims(2, 1)
{ }

/*
 * Constructors for two-dimensional M x N matrices.
 */
template <typename T, typename Syn>
matrix<T,Syn>::matrix(unsigned long M, unsigned long N) 
 : array<T,Syn>(M*N),
   _dims(2)
{
   _dims[0] = M;
   _dims[1] = N;
}

template <typename T, typename Syn>
matrix<T,Syn>::matrix(unsigned long M, unsigned long N, const T& t)
 : array<T,Syn>(M*N, t),
   _dims(2)
{
   _dims[0] = M;
   _dims[1] = N;
}

/*
 * Constructors for multi-dimensional matrices.
 */
template <typename T, typename Syn>
matrix<T,Syn>::matrix(const array<unsigned long>& dims)
 : array<T,Syn>(matrix<T,Syn>::dims_matrix_size(dims)),
   _dims(dims)
{ }
   
template <typename T, typename Syn>
matrix<T,Syn>::matrix(const array<unsigned long>& dims, const T& t)
 : array<T,Syn>(matrix<T,Syn>::dims_matrix_size(dims), t),
   _dims(dims)
{ }

/*
 * Copy constructors.
 */
template <typename T, typename Syn>
matrix<T,Syn>::matrix(const matrix<T,Syn>& m)
 : array<T,Syn>(),
   _dims()
{
   auto_read_lock<const Syn> rlock(m);
   /* copy array */
   matrix<T,Syn> m_copy(m._size, 1);
   const T* m_data = m._data;
   for (unsigned long n = 0; n < m._size; n++)
      m_copy._data[n] = m_data[n];
   /* copy dimensions */
   _dims = m._dims;
   /* take ownership of array copy */
   array<T,Syn>::swap(*this, m_copy);
}

template <typename T, typename Syn>
template <typename U, typename S>
matrix<T,Syn>::matrix(const matrix<U,S>& m)
 : array<T,Syn>(),
   _dims()
{
   auto_read_lock<const S> rlock(m);
   /* copy array */
   matrix<T,Syn> m_copy(m._size, 1);
   const U* m_data = m._data;
   for (unsigned long n = 0; n < m._size; n++)
      m_copy._data[n] = static_cast<T>(m_data[n]);
   /* copy dimensions */
   _dims = m._dims;
   /* take ownership of array copy */
   array<T,Syn>::swap(*this, m_copy);
}

/*
 * Destructor.
 */
template <typename T, typename Syn>
matrix<T,Syn>::~matrix() {
   /* do nothing */
}
   
/***************************************************************************
 * Named constructors.
 ***************************************************************************/
 
/*
 * Empty matrix (0 x 0).
 */
template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::empty() {
   return matrix<T,Syn>();
}

/*
 * Identity matrix.
 * Create and return an N x N identity matrix.
 */
template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::eye(unsigned long N) {
   matrix<T,Syn> m(N, N);
   const T t_one(1);
   for (unsigned long n = 0; n < m._size; n += (N+1))
      m._data[n] = t_one;
   return m;
}

/*
 * Zeros matrix.
 * Create and return an M x N matrix of zeros.
 */
template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::zeros(unsigned long M, unsigned long N) {
   return matrix<T,Syn>(M, N);
}

/*
 * Zeros matrix.
 * Create and return a matrix of zeros with the given dimensions.
 */
template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::zeros(const array<unsigned long>& dims) {
   return matrix<T,Syn>(dims);
}

/*
 * Ones matrix.
 * Create and return an M x N matrix of ones.
 */
template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::ones(unsigned long M, unsigned long N) {
   return matrix<T,Syn>(M, N, T(1));
}

/*
 * Ones matrix.
 * Create and return a matrix of ones with the given dimensions.
 */
template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::ones(const array<unsigned long>& dims) {
   return matrix<T,Syn>(dims, T(1));
}

/*
 * Random matrix.
 * Create and return an M x N matrix of random numbers drawn from the
 * given random number generator.  If no generator is specified, then
 * the matrix elements are uniformly distributed random numbers in the 
 * interval [0, 1] (inclusive).
 */ 
template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::random(
   unsigned long M, unsigned long N)
{
   rand_gen_uniform<T> r;
   return matrix<T,Syn>::random(M, N, r);
}

template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::random(
   unsigned long M, unsigned long N, rand_gen<T>& r)
{
   matrix<T,Syn> m(M, N);
   for (unsigned long n = 0; n < m._size; n++)
      m._data[n] = r.generate();
   return m;
}

/*
 * Random matrix.
 * Create and return a matrix with the given dimensions in which each 
 * element is a random number drawn from the given random number generator.
 * If no generator is specified, then the matrix elements are uniformly 
 * distributed random numbers in the interval [0, 1] (inclusive).
 */
template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::random(
   const array<unsigned long>& dims)
{
   rand_gen_uniform<T> r;
   return matrix<T,Syn>::random(dims, r);
}

template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::random(
   const array<unsigned long>& dims, rand_gen<T>& r)
{
   matrix<T,Syn> m(dims);
   for (unsigned long n = 0; n < m._size; n++)
      m._data[n] = r.generate();
   return m;
}

/*
 * Ramp.
 * Create and return a vector containing values from start to end at an 
 * increment of step.  The default step is one.
 */
template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::ramp(const T& start, const T& end) {
   const T step(1);
   return matrix<T,Syn>::ramp(start, step, end);
}

template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::ramp(const T& start, const T& step, const T& end) {
   /* check for nonzero step */
   if (step == T())
      throw ex_invalid_argument(
         "step must be nonzero"
      );
   /* compute number of vector elements */
   const T space = (end - start) / step;
   unsigned long n_els =
      (space >= T()) ? (static_cast<unsigned long>(space) + 1) : 0;
   /* create vector */
   matrix<T,Syn> m(n_els, 1);
   T curr = start;
   for (unsigned long n = 0; n < n_els; n++) {
      m._data[n] = curr;
      curr += step;
   }
   return m;
}

/*
 * Diagonal matrix.
 * Create and return a d-dimensional square matrix given its diagonal.
 */
template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::diagonal(
   const matrix<T,Syn>& dgnl, unsigned long d)
{
   auto_read_lock<const Syn> rlock(dgnl);
   /* create dimensions and allocate matrix */
   array<unsigned long> dims(d, dgnl._size);
   matrix<T,Syn> m(dims);
   /* compute increment size along diagonal */
   array<unsigned long> sizes = matrix<T,Syn>::dims_sizes(dims);
   unsigned long step = 0;
   for (unsigned long n = 0; n < d; n++)
      step += sizes[n];
   /* set diagonal elements */
   for (unsigned long n = 0, pos = 0; n < m._size; n++, pos += step)
      m._data[pos] = dgnl._data[n];
   return m;
}

/***************************************************************************
 * Copy.
 ***************************************************************************/

/*
 * Copy the contents of a matrix.
 * Note: matrix is not locked by this function.
 */
template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::copy(const matrix<T,Syn>& m) {
   matrix<T,Syn> m_copy(m._dims);
   const T* m_data = m._data;
   for (unsigned long n = 0; n < m._size; n++)
      m_copy._data[n] = m_data[n];
   return m_copy;
}

/*
 * Copy the portion of the source matrix that overlaps the destination
 * matrix into the destination matrix.
 * Note: matrices are not locked by this function.
 */
template <typename T, typename Syn>
void matrix<T,Syn>::copy_overlap(
   const matrix<T,Syn>& src, matrix<T,Syn>& dst)
{
   if ((src._size > 0) && (dst._size > 0)) {
      /* compute start and end of block shared between matrices */
      unsigned long n_dims_src = src._dims.size();
      unsigned long n_dims_dst = dst._dims.size();
      unsigned long n_dims =
         (n_dims_src > n_dims_dst) ? n_dims_src : n_dims_dst;
      array<unsigned long> dims_src =
         matrix<T,Syn>::extend_dims(src._dims, n_dims);
      array<unsigned long> dims_dst =
         matrix<T,Syn>::extend_dims(dst._dims, n_dims);
      array<unsigned long> start(n_dims);
      array<unsigned long> end(n_dims);
      for (unsigned long n = 0; n < n_dims; n++) {
         end[n] =
            ((dims_src[n] < dims_dst[n]) ? dims_src[n] : dims_dst[n]) - 1;
      }
      /* compute indices of block in original and resized matrix */
      array<unsigned long> inds_src = matrix<T,Syn>::linear_indices(
         dims_src, start, end
      );
      array<unsigned long> inds_dst = matrix<T,Syn>::linear_indices(
         dims_dst, start, end
      );
      /* copy shared block */
      unsigned long n_inds = inds_src.size();
      for (unsigned long n = 0; n < n_inds; n++)
         dst._data[inds_dst[n]] = src._data[inds_src[n]];
   }
}

/***************************************************************************
 * Swap.
 ***************************************************************************/

/*
 * Swap the contents of two matrices.
 * Note: matrices are not locked by this function.
 */
template <typename T, typename Syn>
void matrix<T,Syn>::swap(matrix<T,Syn>& m0, matrix<T,Syn>& m1) {
   /* swap dimensions */
   array<unsigned long> temp_dims = m0._dims;
   m0._dims = m1._dims;
   m1._dims = temp_dims;
   /* swap data */
   array<T,Syn>::swap(m0, m1);
}

/***************************************************************************
 * Conversion to matrix.
 ***************************************************************************/

/*
 * Convert an array<T,Syn> object into a matrix<T,Syn> object.
 * The original array<T,Syn> object is destroyed in the process.
 */
template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::to_matrix(array<T,Syn>& a) {
   matrix<T,Syn> m;
   array<T,Syn>::swap(m, a);
   m._dims[0] = m._size;
   m._dims[1] = 1;
   return m;
}

/***************************************************************************
 * Serialization.
 ***************************************************************************/

/*
 * Serialize.
 */
template <typename T, typename Syn>
void matrix<T,Syn>::serialize(
   serial_output_stream& s, const serializer<T>& slzr) const
{
   auto_read_lock<const Syn> rlock(*this);
   _dims.serialize(s);
   this->array<T,Syn>::serialize(s, slzr);
}

/*
 * Deserialize.
 */
template <typename T, typename Syn>
auto_ptr< matrix<T,Syn> > matrix<T,Syn>::deserialize(
   serial_input_stream& s, const serializer<T>& slzr)
{
   /* deserialize dimensions and data */
   auto_ptr< array<unsigned long> > dims = array<unsigned long>::deserialize(s);
   auto_ptr< array<T,Syn> > data = array<T,Syn>::deserialize(s, slzr);
   /* allocate matrix and swap in dimensions, data */
   auto_ptr< matrix<T,Syn> > m(new matrix<T,Syn>());
   array<unsigned long>::swap(m->_dims, *dims);
   array<T,Syn>::swap(*m, *data);
   return m;
}

/***************************************************************************
 * Element reference.
 ***************************************************************************/

/*
 * Element reference (using two indices).
 */
template <typename T, typename Syn>
T& matrix<T,Syn>::operator()(unsigned long x, unsigned long y) {
   auto_read_lock<const Syn> rlock(*this);
   /* compute linear index */
   unsigned long n = matrix<T,Syn>::linear_index(_dims, x, y);
   #if MATH__MATRICES__MATRIX__CHECK_BOUNDS
      /* perform bounds check */
      if (n >= this->_size)
         throw ex_index_out_of_bounds(
            "matrix index out of bounds", 
            n
         );
   #endif
   /* retrieve element */
   return this->_data[n];
}

template <typename T, typename Syn>
const T& matrix<T,Syn>::operator()(unsigned long x, unsigned long y) const {
   auto_read_lock<const Syn> rlock(*this);
   /* compute linear index */
   unsigned long n = matrix<T,Syn>::linear_index(_dims, x, y);
   #if MATH__MATRICES__MATRIX__CHECK_BOUNDS
      /* perform bounds check */
      if (n >= this->_size)
         throw ex_index_out_of_bounds(
            "matrix index out of bounds", 
            n
         );
   #endif
   /* retrieve element */
   return this->_data[n];
}
 
/*
 * Element reference (using multiple indices).
 */
template <typename T, typename Syn>
T& matrix<T,Syn>::operator()(const array<unsigned long>& i) {
   auto_read_lock<const Syn> rlock(*this);
   /* compute linear index */
   unsigned long n = matrix<T,Syn>::linear_index(_dims, i);
   #if MATH__MATRICES__MATRIX__CHECK_BOUNDS
      /* perform bounds check */
      if (n >= this->_size)
         throw ex_index_out_of_bounds(
            "matrix index out of bounds", 
            n
         );
   #endif
   /* retrieve element */
   return this->_data[n];
}

template <typename T, typename Syn>
const T& matrix<T,Syn>::operator()(const array<unsigned long>& i) const {
   auto_read_lock<const Syn> rlock(*this);
   /* compute linear index */
   unsigned long n = matrix<T,Syn>::linear_index(_dims, i);
   #if MATH__MATRICES__MATRIX__CHECK_BOUNDS
      /* perform bounds check */
      if (n >= this->_size)
         throw ex_index_out_of_bounds(
            "matrix index out of bounds", 
            n
         );
   #endif
   /* retrieve element */
   return this->_data[n];
}

/***************************************************************************
 * Direct element access.
 ***************************************************************************/

/*
 * Direct element access.
 * Return a pointer to the memory containing the matrix elements.
 * This should be used with extreme caution.
 */
template <typename T, typename Syn>
T* matrix<T,Syn>::data() {
   auto_read_lock<const Syn> rlock(*this);
   return this->_data;
}

template <typename T, typename Syn>
const T* matrix<T,Syn>::data() const {
   auto_read_lock<const Syn> rlock(*this);
   return this->_data;
}

/***************************************************************************
 * Subarray.
 ***************************************************************************/

/*
 * Subarray.
 * Return vector of elements in the specified linear index range.
 */
template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::subarray(
   unsigned long start,
   unsigned long end) const
{
   array<T,Syn> a = this->array<T,Syn>::subarray(start, end);
   return matrix<T,Syn>::to_matrix(a);
}

template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::subarray(
   unsigned long start, 
   unsigned long step,
   unsigned long end) const
{
   array<T,Syn> a = this->array<T,Syn>::subarray(start, step, end);
   return matrix<T,Syn>::to_matrix(a);
}

template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::subarray(
   const array<unsigned long>& inds) const
{
   array<T,Syn> a = this->array<T,Syn>::subarray(inds);
   return matrix<T,Syn>::to_matrix(a);
}

/***************************************************************************
 * Submatrix.
 ***************************************************************************/
 
/*
 * Submatrix.
 * Return matrix of elements in the specified range.
 */
template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::submatrix(
   const array<unsigned long>& start,
   const array<unsigned long>& end) const
{
   array< array<unsigned long> > range = matrix<T,Syn>::range_indices(
      start, end
   );
   return this->submatrix(range);
}

template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::submatrix(
   const array<unsigned long>& start,
   const array<unsigned long>& step,
   const array<unsigned long>& end) const
{
   array< array<unsigned long> > range = matrix<T,Syn>::range_indices(
      start, step, end
   );
   return this->submatrix(range);
}

template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::submatrix(
   const array< array<unsigned long> >& range) const
{
   auto_read_lock<const Syn> rlock(*this);
   /* check arguments */
   unsigned long n_dims = _dims.size();
   if (range.size() != n_dims)
      throw ex_invalid_argument(
         "specified range for submatrix has incorrent number of dimensions"
      );
   /* compute dimensions of resulting matrix */
   array<unsigned long> subm_dims(n_dims);
   for (unsigned long n = 0; n < n_dims; n++)
      subm_dims[n] = range[n].size();
   /* translate range to linear indices */
   array<unsigned long> inds = matrix<T,Syn>::linear_indices(_dims, range);
   /* grab desired elements */
   matrix<T,Syn> m(subm_dims);
   for (unsigned long n = 0; n < m._size; n++) {
      unsigned long ind = inds[n];
      #if MATH__MATRICES__MATRIX__CHECK_BOUNDS
         /* perform bounds check */
         if (ind >= this->_size)
            throw ex_index_out_of_bounds(
               "matrix index out of bounds", 
               ind
            );
      #endif
      m._data[n] = this->_data[ind];
   }
   return m;
}

/***************************************************************************
 * Subassign.
 ***************************************************************************/

/*
 * Subarray assignment.
 * Assign value(s) to elements in the specified linear index range.
 */
template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::subassign(
   unsigned long start,
   unsigned long end,
   const T& t)
{
   return this->subassign(start, 1, end, t);
}

template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::subassign(
   unsigned long start,
   unsigned long end,
   const matrix<T,Syn>& values)
{
   return this->subassign(start, 1, end, values);
}

template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::subassign(
   unsigned long start,
   unsigned long step,
   unsigned long end,
   const T& t)
{
   auto_write_lock<const Syn> wlock(*this);
   /* check step size */
   if (step == 0)
      throw ex_invalid_argument("step must be nonzero");
   /* check start and end indices */
   if (start >= (this->_size))
      throw ex_index_out_of_bounds("start index exceeds matrix size", start);
   if (end >= (this->_size))
      throw ex_index_out_of_bounds("end index exceeds matrix size", end);
   /* assign value to elements */
   const T t_copy(t);
   for (unsigned long n = start; n <= end; n += step)
      this->_data[n] = t_copy;
   return *this;
}

template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::subassign(
   unsigned long start,
   unsigned long step,
   unsigned long end,
   const matrix<T,Syn>& values)
{
   auto_read_write_lock<const Syn> rwlock(values, *this);
   /* check step size */
   if (step == 0)
      throw ex_invalid_argument("step must be nonzero");
   /* check start and end indices */
   if (start >= (this->_size))
      throw ex_index_out_of_bounds("start index exceeds matrix size", start);
   if (end >= (this->_size))
      throw ex_index_out_of_bounds("end index exceeds matrix size", end);
   /* check number of indices */
   unsigned long n_inds = (start <= end) ? ((end - start)/step + 1) : 0;
   if (n_inds != values._size)
      throw ex_invalid_argument(
         "incorrect number of values given for specified range"
      );
   /* assign values to elements */
   for (unsigned long n = start, n_v = 0; n <= end; n += step, n_v++)
      this->_data[n] = values._data[n_v];
   return *this;
}

template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::subassign(
   const array<unsigned long>& inds,
   const T& t)
{
   auto_write_lock<const Syn> wlock(*this);
   const T t_copy(t);
   unsigned long n_inds = inds.size();
   for (unsigned long n = 0; n < n_inds; n++) {
      unsigned long ind = inds[n];
      #if MATH__MATRICES__MATRIX__CHECK_BOUNDS
         /* perform bounds check */
         if (ind >= this->_size)
            throw ex_index_out_of_bounds(
               "matrix index out of bounds", 
               ind
            );
      #endif
      this->_data[ind] = t_copy;
   }
   return *this;
}

template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::subassign(
   const array<unsigned long>& inds,
   const matrix<T,Syn>& values)
{
   auto_read_write_lock<const Syn> rwlock(values, *this);
   /* check number of indices */
   unsigned long n_inds = inds.size();
   if (n_inds != values._size)
      throw ex_invalid_argument(
         "incorrect number of values given for specified range"
      );
   /* assign values to elements */
   for (unsigned long n = 0; n < n_inds; n++) {
      unsigned long ind = inds[n];
      #if MATH__MATRICES__MATRIX__CHECK_BOUNDS
         /* perform bounds check */
         if (ind >= this->_size)
            throw ex_index_out_of_bounds(
               "matrix index out of bounds", 
               ind
            );
      #endif
      this->_data[ind] = values._data[n];
   }
   return *this;
}

/*
 * Submatrix assignment.
 * Assign value(s) to elements in the specified range.
 */
template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::subassign(
   const array<unsigned long>& start,
   const array<unsigned long>& end,
   const T& t)
{
   array< array<unsigned long> > range = matrix<T,Syn>::range_indices(
      start, end
   );
   return this->subassign(range, t);
}

template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::subassign(
   const array<unsigned long>& start,
   const array<unsigned long>& end,
   const matrix<T,Syn>& values)
{
   array< array<unsigned long> > range = matrix<T,Syn>::range_indices(
      start, end
   );
   return this->subassign(range, values);
}

template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::subassign(
   const array<unsigned long>& start,
   const array<unsigned long>& step,
   const array<unsigned long>& end,
   const T& t)
{
   array< array<unsigned long> > range = matrix<T,Syn>::range_indices(
      start, step, end
   );
   return this->subassign(range, t);
}
  
template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::subassign(
   const array<unsigned long>& start,
   const array<unsigned long>& step,
   const array<unsigned long>& end,
   const matrix<T,Syn>& values)
{
   array< array<unsigned long> > range = matrix<T,Syn>::range_indices(
      start, step, end
   );
   return this->subassign(range, values);
}
 
template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::subassign(
   const array< array<unsigned long> >& range,
   const T& t)
{
   auto_write_lock<const Syn> wlock(*this);
   const T t_copy(t);
   /* translate range to linear indices */
   array<unsigned long> inds = matrix<T,Syn>::linear_indices(_dims, range);
   /* assign values to elements */
   unsigned long n_inds = inds.size();
   for (unsigned long n = 0; n < n_inds; n++) {
      unsigned long ind = inds[n];
      #if MATH__MATRICES__MATRIX__CHECK_BOUNDS
         /* perform bounds check */
         if (ind >= this->_size)
            throw ex_index_out_of_bounds(
               "matrix index out of bounds", 
               ind
            );
      #endif
      this->_data[ind] = t_copy;
   }
   return *this;
}

template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::subassign(
   const array< array<unsigned long> >& range,
   const matrix<T,Syn>& values)
{
   auto_read_write_lock<const Syn> rwlock(values, *this);
   /* translate range to linear indices */
   array<unsigned long> inds = matrix<T,Syn>::linear_indices(_dims, range);
   /* check number of indices */
   unsigned long n_inds = inds.size();
   if (n_inds != values._size)
      throw ex_invalid_argument(
         "incorrect number of values given for specified range"
      );
   /* assign values to elements */
   for (unsigned long n = 0; n < n_inds; n++) {
      unsigned long ind = inds[n];
      #if MATH__MATRICES__MATRIX__CHECK_BOUNDS
         /* perform bounds check */
         if (ind >= this->_size)
            throw ex_index_out_of_bounds(
               "matrix index out of bounds", 
               ind
            );
      #endif
      this->_data[ind] = values._data[n];
   }
   return *this;
}

/***************************************************************************
 * Fill.
 ***************************************************************************/

/*
 * Fill with value.
 */
template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::fill(const T& t) {
   const T t_copy(t);
   auto_write_lock<const Syn> wlock(*this);
   for (unsigned long n = 0; n < this->_size; n++)
      this->_data[n] = t_copy;
   return *this;
}

/***************************************************************************
 * I/O.
 ***************************************************************************/

/*
 * Formatted output to stream.
 */
template <typename T, typename Syn>
ostream& operator<<(ostream& os, const matrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   /* output matrix dimensions */
   unsigned long n_dims = m._dims.size();
   m.output_classname(os);
   os << ": ";
   if (n_dims == 0) {
      os << "[]";
   } else if (n_dims == 1) {
      os << "[" << m._dims[0] << "]";
   } else {
      os << "[";
      for (unsigned long n = 0; n < (n_dims-1); n++)
         os << m._dims[n] << " x ";
      os << m._dims[n_dims-1] << "]";
   }
   os << " (" << m._size << " element(s))";
   if (m._size > 0)
      os << "\n";
   /* initialize indices counter */ 
   array<unsigned long> indices(n_dims);
   for (unsigned long n = 0; n < n_dims; n++)
      indices[n] = 0;
   /* output matrix elements */
   unsigned long ind = 0;
   for (unsigned long n = 0; n < m._size; n++) {
      /* output level break if needed */
      if ((ind+2) < n_dims) {
         os << '(';
         for (ind = 0; ind < (n_dims-2); ind++)
            os << indices[ind] << ',';
         os << ":,:) =\n";
      }
      /* output leading spaces if needed */
      ind = n_dims - 1;
      if (indices[ind] == 0)
         os << "   ";
      /* output element */
      m.output_element(os, n);
      /* increment indices counter */
      indices[ind]++;
      /* output line break if needed */
      if ((n+1) < m._size) {
         if (indices[ind] == m._dims[ind])
            os << "\n";
         else
            os << " ";
      }
      /* update indices counter */
      while ((indices[ind] == m._dims[ind]) && (ind > 0)) {
         indices[ind] = 0;
         ind--;
         indices[ind]++;
      }
   }
   return os;
}

/*
 * Output classname to stream.
 */
template <typename T, typename Syn>
ostream& matrix<T,Syn>::output_classname(ostream& os) const {
   os << "matrix";
   return os;
}

/*
 * Output matrix element at given index to stream.
 */
template <typename T, typename Syn>
ostream& matrix<T,Syn>::output_element(ostream& os, unsigned long n) const {
   /* output to string */
   ostringstream s;
   s << io::streams::iomanip::setw(10)
     << io::streams::iomanip::setprecision(4)
     << io::streams::iomanip::setiosflags(ios::right)
     << this->_data[n];
   /* output string */
   os << s.str();
   return os;
}

/***************************************************************************
 * Assignment operators: matrix-real.
 ***************************************************************************/

template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::operator=(const T& t) {
   matrix<T,Syn> m(t);
   auto_write_lock<const Syn> wlock(*this);
   matrix<T,Syn>::swap(*this, m);
   return *this;
}

template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::operator+=(const T& t) {
   const T t_copy(t);
   auto_write_lock<const Syn> wlock(*this);
   for (unsigned long n = 0; n < this->_size; n++)
      this->_data[n] += t_copy;
   return *this;
}

template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::operator-=(const T& t) {
   const T t_copy(t);
   auto_write_lock<const Syn> wlock(*this);
   for (unsigned long n = 0; n < this->_size; n++)
      this->_data[n] -= t_copy;
   return *this;
}

template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::operator*=(const T& t) {
   const T t_copy(t);
   auto_write_lock<const Syn> wlock(*this);
   for (unsigned long n = 0; n < this->_size; n++)
      this->_data[n] *= t_copy;
   return *this;
}

template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::operator/=(const T& t) {
   const T t_copy(t);
   auto_write_lock<const Syn> wlock(*this);
   for (unsigned long n = 0; n < this->_size; n++)
      this->_data[n] /= t_copy;
   return *this;
}

/***************************************************************************
 * Assignment operators: matrix-matrix.
 ***************************************************************************/

template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::operator=(const matrix<T,Syn>& m) {
   matrix<T,Syn> m_copy(m);
   auto_write_lock<const Syn> wlock(*this);
   matrix<T,Syn>::swap(*this, m_copy);
   return *this;
}

template <typename T, typename Syn>
template <typename U, typename S> 
matrix<T,Syn>& matrix<T,Syn>::operator=(const matrix<U,S>& m) {
   matrix<T,Syn> m_copy(m);
   auto_write_lock<const Syn> wlock(*this);
   matrix<T,Syn>::swap(*this, m_copy);
   return *this;
}

template <typename T, typename Syn>
template <typename U, typename S> 
matrix<T,Syn>& matrix<T,Syn>::operator+=(const matrix<U,S>& m) {
   auto_read_write_lock<const synchronizable> rwlock(m, *this);
   matrix<T,Syn>::assert_dims_equal(_dims, m._dims);
   for (unsigned long n = 0; n < this->_size; n++)
      this->_data[n] += m._data[n];
   return *this;
}

template <typename T, typename Syn>
template <typename U, typename S> 
matrix<T,Syn>& matrix<T,Syn>::operator-=(const matrix<U,S>& m) {
   auto_read_write_lock<const synchronizable> rwlock(m, *this);
   matrix<T,Syn>::assert_dims_equal(_dims, m._dims);
   for (unsigned long n = 0; n < this->_size; n++)
      this->_data[n] -= m._data[n];
   return *this;
}

template <typename T, typename Syn>
template <typename U, typename S> 
matrix<T,Syn>& matrix<T,Syn>::operator*=(const matrix<U,S>& m) {
   auto_read_write_lock<const synchronizable> rwlock(m, *this);
   /* compute result dimensions */ 
   array<unsigned long> dims = matrix<T,Syn>::dims_mult(_dims, m._dims);
   unsigned long n_dims_left       = _dims.size();
   unsigned long n_dims_right      = m._dims.size();
   unsigned long n_dims_from_left  = (n_dims_left > 2)  ? (n_dims_left - 1)  : 1;
   unsigned long n_dims_from_right = (n_dims_right > 2) ? (n_dims_right - 1) : 1;
   unsigned long n_dims_result     = n_dims_from_left + n_dims_from_right;
   /* compute # rows in left argument */
   unsigned long rowsL = 1;
   for (unsigned long n = 0; n < n_dims_from_left; n++)
      rowsL *= dims[n];
   /* compute # cols in right argument */
   unsigned long colsR = 1;
   for (unsigned long n = n_dims_from_left; n < n_dims_result; n++)
      colsR *= dims[n];
   /* compute shared dimension */
   unsigned long colsL_rowsR = (n_dims_right == 0) ? 0 : m._dims[0];
   /* multiply */
   matrix<T,Syn> m_result(dims);
   T*       M  = m_result._data;
   const T* ML = this->_data;
   const U* MR = m._data;
   const U* MR_ptr;
   for (unsigned long r = 0; r < rowsL; r++) {
      M  += colsR;
      ML += colsL_rowsR;
      MR_ptr = MR;
      for (const T* ML_ptr = (ML - colsL_rowsR); ML_ptr < ML; ML_ptr++) {
         for (T* M_ptr = (M - colsR); M_ptr < M; M_ptr++) {
            (*M_ptr) += (*ML_ptr) * (*MR_ptr);
            MR_ptr++;
         }
      }
   }
   /* take ownership of result */
   matrix<T,Syn>::swap(*this, m_result);
   return *this;
}

template <typename T, typename Syn>
template <typename U, typename S> 
matrix<T,Syn>& matrix<T,Syn>::operator/=(const matrix<U,S>& m) {
   (*this) *= inv(m);
   return *this;
} 

/***************************************************************************
 * Binary matrix-real operators.
 ***************************************************************************/
 
template <typename T, typename Syn>
matrix<T,Syn> operator+(const matrix<T,Syn>& m, const T& t) {
   return (matrix<T,Syn>(m) += t);
}

template <typename T, typename Syn>
matrix<T,Syn> operator-(const matrix<T,Syn>& m, const T& t) {
   return (matrix<T,Syn>(m) -= t);
}

template <typename T, typename Syn>
matrix<T,Syn> operator*(const matrix<T,Syn>& m, const T& t) {
   return (matrix<T,Syn>(m) *= t);
}

template <typename T, typename Syn>
matrix<T,Syn> operator/(const matrix<T,Syn>& m, const T& t) {
   return (matrix<T,Syn>(m) /= t);
}

template <typename T, typename Syn>
matrix<T,Syn> operator+(const T& t, const matrix<T,Syn>& m) {
   return (matrix<T,Syn>(m) += t);
}

template <typename T, typename Syn>
matrix<T,Syn> operator-(const T& t, const matrix<T,Syn>& m) {
   matrix<T,Syn> m_result = (-m);
   m_result += t;
   return m_result;
}

template <typename T, typename Syn>
matrix<T,Syn> operator*(const T& t, const matrix<T,Syn>& m) {
   return (matrix<T,Syn>(m) *= t);
}

template <typename T, typename Syn>
matrix<T,Syn> operator/(const T& t, const matrix<T,Syn>& m) {
   matrix<T,Syn> m_result = inv(m);
   m_result *= t;
   return m_result;
}

/***************************************************************************
 * Binary matrix-matrix operators.
 ***************************************************************************/

template <typename T, typename Syn>
matrix<T,Syn> operator+(const matrix<T,Syn>& m0, const matrix<T,Syn>& m1) {
   return (matrix<T,Syn>(m0) += m1);
}

template <typename T, typename Syn>
matrix<T,Syn> operator-(const matrix<T,Syn>& m0, const matrix<T,Syn>& m1) {
   return (matrix<T,Syn>(m0) -= m1);
}

template <typename T, typename Syn>
matrix<T,Syn> operator*(const matrix<T,Syn>& m0, const matrix<T,Syn>& m1) {
   return (matrix<T,Syn>(m0) *= m1);
}

template <typename T, typename Syn>
matrix<T,Syn> operator/(const matrix<T,Syn>& m0, const matrix<T,Syn>& m1) {
   return (matrix<T,Syn>(m0) /= m1);
}

/***************************************************************************
 * Binary matrix-real comparators.
 ***************************************************************************/

template <typename T, typename Syn>
matrix<bool,Syn> operator==(const matrix<T,Syn>& m, const T& t) {
   const T t_copy(t);
   auto_read_lock<const Syn> rlock(m);
   matrix<bool,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result[n] = (m._data[n] == t_copy);
   return m_result;
}

template <typename T, typename Syn>
matrix<bool,Syn> operator!=(const matrix<T,Syn>& m, const T& t) {
   const T t_copy(t);
   auto_read_lock<const Syn> rlock(m);
   matrix<bool,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result[n] = (m._data[n] != t_copy);
   return m_result;
}

template <typename T, typename Syn>
matrix<bool,Syn> operator<(const matrix<T,Syn>& m, const T& t) {
   const T t_copy(t);
   auto_read_lock<const Syn> rlock(m);
   matrix<bool,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result[n] = (m._data[n] < t_copy);
   return m_result;
}

template <typename T, typename Syn>
matrix<bool,Syn> operator>(const matrix<T,Syn>& m, const T& t) {
   const T t_copy(t);
   auto_read_lock<const Syn> rlock(m);
   matrix<bool,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result[n] = (m._data[n] > t_copy);
   return m_result;
}

template <typename T, typename Syn>
matrix<bool,Syn> operator<=(const matrix<T,Syn>& m, const T& t) {
   const T t_copy(t);
   auto_read_lock<const Syn> rlock(m);
   matrix<bool,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result[n] = (m._data[n] <= t_copy);
   return m_result;
}

template <typename T, typename Syn>
matrix<bool,Syn> operator>=(const matrix<T,Syn>& m, const T& t) {
   const T t_copy(t);
   auto_read_lock<const Syn> rlock(m);
   matrix<bool,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result[n] = (m._data[n] >= t_copy);
   return m_result;
}

template <typename T, typename Syn>
matrix<bool,Syn> operator==(const T& t, const matrix<T,Syn>& m) {
   return (m == t);
}

template <typename T, typename Syn>
matrix<bool,Syn> operator!=(const T& t, const matrix<T,Syn>& m) {
   return (m != t);
}

template <typename T, typename Syn>
matrix<bool,Syn> operator<(const T& t, const matrix<T,Syn>& m) {
   return (m > t);
}

template <typename T, typename Syn>
matrix<bool,Syn> operator>(const T& t, const matrix<T,Syn>& m) {
   return (m < t);
}

template <typename T, typename Syn>
matrix<bool,Syn> operator<=(const T& t, const matrix<T,Syn>& m) {
   return (m >= t);
}

template <typename T, typename Syn>
matrix<bool,Syn> operator>=(const T& t, const matrix<T,Syn>& m) {
   return (m <= t);
}

/***************************************************************************
 * Binary matrix-complex comparators.
 ***************************************************************************/

template <typename T, typename Syn>
matrix<bool,Syn> operator==(const matrix<T,Syn>& m, const complex<T>& z) {
   const complex<T> z_copy(z);
   auto_read_lock<const Syn> rlock(m);
   matrix<bool,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result[n] = (m._data[n] == z_copy);
   return m_result;
}

template <typename T, typename Syn>
matrix<bool,Syn> operator!=(const matrix<T,Syn>& m, const complex<T>& z) {
   const complex<T> z_copy(z);
   auto_read_lock<const Syn> rlock(m);
   matrix<bool,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result[n] = (m._data[n] != z_copy);
   return m_result;
}

template <typename T, typename Syn>
matrix<bool,Syn> operator<(const matrix<T,Syn>& m, const complex<T>& z) {
   const complex<T> z_copy(z);
   auto_read_lock<const Syn> rlock(m);
   matrix<bool,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result[n] = (m._data[n] < z_copy);
   return m_result;
}

template <typename T, typename Syn>
matrix<bool,Syn> operator>(const matrix<T,Syn>& m, const complex<T>& z) {
   const complex<T> z_copy(z);
   auto_read_lock<const Syn> rlock(m);
   matrix<bool,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result[n] = (m._data[n] > z_copy);
   return m_result;
}

template <typename T, typename Syn>
matrix<bool,Syn> operator<=(const matrix<T,Syn>& m, const complex<T>& z) {
   const complex<T> z_copy(z);
   auto_read_lock<const Syn> rlock(m);
   matrix<bool,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result[n] = (m._data[n] <= z_copy);
   return m_result;
}

template <typename T, typename Syn>
matrix<bool,Syn> operator>=(const matrix<T,Syn>& m, const complex<T>& z) {
   const complex<T> z_copy(z);
   auto_read_lock<const Syn> rlock(m);
   matrix<bool,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result[n] = (m._data[n] >= z_copy);
   return m_result;
}

template <typename T, typename Syn>
matrix<bool,Syn> operator==(const complex<T>& z, const matrix<T,Syn>& m) {
   return (m == z);
}

template <typename T, typename Syn>
matrix<bool,Syn> operator!=(const complex<T>& z, const matrix<T,Syn>& m) {
   return (m != z);
}

template <typename T, typename Syn>
matrix<bool,Syn> operator<(const complex<T>& z, const matrix<T,Syn>& m) {
   return (m > z);
}

template <typename T, typename Syn>
matrix<bool,Syn> operator>(const complex<T>& z, const matrix<T,Syn>& m) {
   return (m < z);
}

template <typename T, typename Syn>
matrix<bool,Syn> operator<=(const complex<T>& z, const matrix<T,Syn>& m) {
   return (m >= z);
}

template <typename T, typename Syn>
matrix<bool,Syn> operator>=(const complex<T>& z, const matrix<T,Syn>& m) {
   return (m <= z);
}

/***************************************************************************
 * Binary matrix-matrix comparators.
 ***************************************************************************/

template <typename T, typename Syn>
matrix<bool,Syn> operator==(const matrix<T,Syn>& m0, const matrix<T,Syn>& m1) {
   auto_read_read_lock<const Syn> rrlock(m0, m1);
   matrix<T,Syn>::assert_dims_equal(m0._dims, m1._dims);
   matrix<bool,Syn> m(m0._dims);
   for (unsigned long n = 0; n < m0._size; n++)
      m[n] = (m0._data[n] == m1._data[n]);
   return m;
}

template <typename T, typename Syn>
matrix<bool,Syn> operator!=(const matrix<T,Syn>& m0, const matrix<T,Syn>& m1) {
   auto_read_read_lock<const Syn> rrlock(m0, m1);
   matrix<T,Syn>::assert_dims_equal(m0._dims, m1._dims);
   matrix<bool,Syn> m(m0._dims);
   for (unsigned long n = 0; n < m0._size; n++)
      m[n] = (m0._data[n] != m1._data[n]);
   return m;
}

template <typename T, typename Syn>
matrix<bool,Syn> operator<(const matrix<T,Syn>& m0, const matrix<T,Syn>& m1) {
   auto_read_read_lock<const Syn> rrlock(m0, m1);
   matrix<T,Syn>::assert_dims_equal(m0._dims, m1._dims);
   matrix<bool,Syn> m(m0._dims);
   for (unsigned long n = 0; n < m0._size; n++)
      m[n] = (m0._data[n] < m1._data[n]);
   return m;
}

template <typename T, typename Syn>
matrix<bool,Syn> operator>(const matrix<T,Syn>& m0, const matrix<T,Syn>& m1) {
   auto_read_read_lock<const Syn> rrlock(m0, m1);
   matrix<T,Syn>::assert_dims_equal(m0._dims, m1._dims);
   matrix<bool,Syn> m(m0._dims);
   for (unsigned long n = 0; n < m0._size; n++)
      m[n] = (m0._data[n] > m1._data[n]);
   return m;
}

template <typename T, typename Syn>
matrix<bool,Syn> operator<=(const matrix<T,Syn>& m0, const matrix<T,Syn>& m1) {
   auto_read_read_lock<const Syn> rrlock(m0, m1);
   matrix<T,Syn>::assert_dims_equal(m0._dims, m1._dims);
   matrix<bool,Syn> m(m0._dims);
   for (unsigned long n = 0; n < m0._size; n++)
      m[n] = (m0._data[n] <= m1._data[n]);
   return m;
}

template <typename T, typename Syn>
matrix<bool,Syn> operator>=(const matrix<T,Syn>& m0, const matrix<T,Syn>& m1) {
   auto_read_read_lock<const Syn> rrlock(m0, m1);
   matrix<T,Syn>::assert_dims_equal(m0._dims, m1._dims);
   matrix<bool,Syn> m(m0._dims);
   for (unsigned long n = 0; n < m0._size; n++)
      m[n] = (m0._data[n] >= m1._data[n]);
   return m;
}

/***************************************************************************
 * Unary operators.
 ***************************************************************************/

/*
 * Unary arithmetic operators.
 */
template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::operator+() const {
   matrix<T,Syn> m(*this);
   for (unsigned long n = 0; n < m._size; n++)
      m._data[n] = +(m._data[n]);
   return m;
}

template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::operator-() const {
   matrix<T,Syn> m(*this);
   for (unsigned long n = 0; n < m._size; n++)
      m._data[n] = -(m._data[n]);
   return m;
}

template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::operator++() {
   auto_write_lock<const Syn> wlock(*this);
   for (unsigned long n = 0; n < this->_size; n++)
      ++(this->_data[n]);
   return *this;
}

template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::operator--() {
   auto_write_lock<const Syn> wlock(*this);
   for (unsigned long n = 0; n < this->_size; n++)
      --(this->_data[n]);
   return *this;
}

template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::operator++(int) {
   auto_write_lock<const Syn> wlock(*this);
   matrix<T,Syn> m(_dims);
   for (unsigned long n = 0; n < this->_size; n++)
      m._data[n] = (this->_data[n])++;
   return m;
}

template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::operator--(int) {
   auto_write_lock<const Syn> wlock(*this);
   matrix<T,Syn> m(_dims);
   for (unsigned long n = 0; n < this->_size; n++)
      m._data[n] = (this->_data[n])--;
   return m;
}

/*
 * Unary element logical inversion operator.
 */
template <typename T, typename Syn>
matrix<bool,Syn> matrix<T,Syn>::operator!() const {
   return (*this == T());
}

/*
 * Unary element multiplicative inversion operator.
 */
template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::operator~() const {
   matrix<T,Syn> m(*this);
   const T t_one(1);
   for (unsigned long n = 0; n < m._size; n++)
      m._data[n] = t_one/m._data[n];
   return m;
}

/***************************************************************************
 * Transcendentals and other functions on elements.
 ***************************************************************************/

template <typename T, typename Syn>
matrix<T,Syn> cos(const matrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   matrix<T,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result._data[n] = cos(m._data[n]);
   return m_result;
}

template <typename T, typename Syn>
matrix<T,Syn> sin(const matrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   matrix<T,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result._data[n] = sin(m._data[n]);
   return m_result;
}

template <typename T, typename Syn>
matrix<T,Syn> tan(const matrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   matrix<T,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result._data[n] = tan(m._data[n]);
   return m_result;
}

template <typename T, typename Syn>
matrix<T,Syn> cosh(const matrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   matrix<T,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result._data[n] = cosh(m._data[n]);
   return m_result;
}

template <typename T, typename Syn>
matrix<T,Syn> sinh(const matrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   matrix<T,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result._data[n] = sinh(m._data[n]);
   return m_result;
}

template <typename T, typename Syn>
matrix<T,Syn> tanh(const matrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   matrix<T,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result._data[n] = tanh(m._data[n]);
   return m_result;
}

template <typename T, typename Syn>
matrix<T,Syn> atan(const matrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   matrix<T,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result._data[n] = atan(m._data[n]);
   return m_result;
}

template <typename T, typename Syn>
matrix<T,Syn> asinh(const matrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   matrix<T,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result._data[n] = asinh(m._data[n]);
   return m_result;
}

template <typename T, typename Syn>
matrix<T,Syn> exp(const matrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   matrix<T,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result._data[n] = exp(m._data[n]);
   return m_result;
}

template <typename T, typename Syn>
matrix<T,Syn> pow(const matrix<T,Syn>& m, int p) {
   auto_read_lock<const Syn> rlock(m);
   matrix<T,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result._data[n] = pow(m._data[n], p);
   return m_result;
}

/*
 * Two-parameter arctangent of the corresponding elements.
 */
template <typename T, typename Syn>
matrix<T,Syn> atan2(const matrix<T,Syn>& my, const matrix<T,Syn>& mx) {
   auto_read_read_lock<const Syn> rrlock(my, mx);
   matrix<T,Syn>::assert_dims_equal(my._dims, mx._dims);
   matrix<T,Syn> m_result(my._dims);
   for (unsigned long n = 0; n < my._size; n++)
      m_result._data[n] = atan2(my._data[n], mx._data[n]);
   return m_result;
}

/*
 * Absolute value of elements.
 */
template <typename T, typename Syn>
matrix<T,Syn> abs(const matrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   matrix<T,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result._data[n] = abs(m._data[n]);
   return m_result;
}

/***************************************************************************
 * Convolution, dot product, and element product.
 ***************************************************************************/

/*
 * Convolution.
 */
template <typename T, typename Syn>
matrix<T,Syn> conv(
   const matrix<T,Syn>& m0, const matrix<T,Syn>& m1)
{
   auto_read_read_lock<const Syn> rrlock(m0, m1);
   return matrix<T,Syn>::compute_conv(m0, m1, false, false);
}

/*
 * Convolution.
 * Crop the result to be no larger than the left input matrix.
 */
template <typename T, typename Syn>
matrix<T,Syn> conv_crop(
   const matrix<T,Syn>& m0, const matrix<T,Syn>& m1)
{
   auto_read_read_lock<const Syn> rrlock(m0, m1);
   return matrix<T,Syn>::compute_conv(m0, m1, true, false);
}

/*
 * Convolution.
 * Crop the result to be no larger than the left input matrix and to
 * include only the central portion for which no zero padding was required.
 */
template <typename T, typename Syn>
matrix<T,Syn> conv_crop_strict(
   const matrix<T,Syn>& m0, const matrix<T,Syn>& m1)
{
   auto_read_read_lock<const Syn> rrlock(m0, m1);
   return matrix<T,Syn>::compute_conv(m0, m1, true, true);
}

/*
 * Compute convolution (for 1D matrices).
 */
template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::compute_conv_1D(
   const matrix<T,Syn>& m0,
   const matrix<T,Syn>& m1,
   bool is_cropped,
   bool is_strict)
{
   /* get size of each matrix */
   unsigned long size0 = m0._size;
   unsigned long size1 = m1._size;
   /* compute dimensions of resulting matrix and also     */ 
   /* compute range of position within full result matrix */
   unsigned long size = 0;
   unsigned long pos_start = 0;
   if (!is_cropped) {
      /* set dimensions for full result matrix (start position is zero) */
      size = ((size0 > 0) && (size1 > 0)) ? (size0 + size1 - 1) : 0;
   } else if (!is_strict) {
      /* set dimensions for result matrix no larger than left input */
      size = ((size0 > 0) && (size1 > 0)) ? (size0) : 0;
      /* set start position for result matrix no larger than left input */
      pos_start = size1/2;
   } else {
      /* set dimensions for central portion of result matrix */
      size = ((size0 >= size1) && (size1 > 0)) ? (size0 - size1 + 1) : 0;
      /* set start position for central portion of result matrix */
      pos_start = (size1 > 0) ? (size1 - 1) : 0;
   }
   /* initialize result */
   array<unsigned long> dims(1);
   dims[0] = size;
   matrix<T,Syn> m(dims);
   /* initialize position in result */
   unsigned long pos = pos_start;
   for (unsigned long n = 0; n < size; n++) {
      /* compute range of offset */
      unsigned long offset_min = ((pos + 1) > size0) ? (pos + 1 - size0) : 0;
      unsigned long offset_max = (pos < size1) ? pos : (size1 - 1);
      /* multiply and add corresponing elements */
      unsigned long ind0 = pos - offset_min;
      unsigned long ind1 = offset_min;
      while (ind1 <= offset_max) {
         /* update result value */
         m._data[n] += m0._data[ind0] * m1._data[ind1];
         /* update linear positions */
         ind0--;
         ind1++;
      }
      /* update position */
      pos++;
   }
   return m;
}

/*
 * Compute convolution (for 2D matrices).
 */
template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::compute_conv_2D(
   const matrix<T,Syn>& m0,
   const matrix<T,Syn>& m1,
   const array<unsigned long>& dims0,
   const array<unsigned long>& dims1,
   bool is_cropped,
   bool is_strict)
{
   /* get size of each matrix */
   unsigned long size0_x = dims0[0];
   unsigned long size0_y = dims0[1];
   unsigned long size1_x = dims1[0];
   unsigned long size1_y = dims1[1];
   /* compute dimensions of resulting matrix and also     */ 
   /* compute range of position within full result matrix */
   unsigned long size_x = 0;
   unsigned long size_y = 0;
   unsigned long pos_start_x = 0;
   unsigned long pos_start_y = 0;
   if (!is_cropped) {
      /* set dimensions for full result matrix (start position is zero) */
      size_x = ((size0_x > 0) && (size1_x > 0)) ? (size0_x + size1_x - 1) : 0;
      size_y = ((size0_y > 0) && (size1_y > 0)) ? (size0_y + size1_y - 1) : 0;
   } else if (!is_strict) {
      /* set dimensions for result matrix no larger than left input */
      size_x = ((size0_x > 0) && (size1_x > 0)) ? (size0_x) : 0;
      size_y = ((size0_y > 0) && (size1_y > 0)) ? (size0_y) : 0;
      /* set start position for result matrix no larger than left input */
      pos_start_x = size1_x/2;
      pos_start_y = size1_y/2;
   } else {
      /* set dimensions for central portion of result matrix */
      size_x =
         ((size0_x >= size1_x) && (size1_x > 0)) ? (size0_x - size1_x + 1) : 0;
      size_y =
         ((size0_y >= size1_y) && (size1_y > 0)) ? (size0_y - size1_y + 1) : 0;
      /* set start position for central portion of result matrix */
      pos_start_x = (size1_x > 0) ? (size1_x - 1) : 0;
      pos_start_y = (size1_y > 0) ? (size1_y - 1) : 0;
   }
   unsigned long pos_bound_y = pos_start_y + size_y;
   /* initialize result */
   matrix<T,Syn> m(size_x, size_y);
   if (m._size == 0)
      return m;
   /* initialize position in result */
   unsigned long pos_x = pos_start_x;
   unsigned long pos_y = pos_start_y;
   /* compute initial range of offset_x */
   unsigned long offset_min_x =
      ((pos_x + 1) > size0_x) ? (pos_x + 1 - size0_x) : 0;
   unsigned long offset_max_x =
      (pos_x < size1_x) ? pos_x : (size1_x - 1);
   unsigned long ind0_start_x = (pos_x - offset_min_x) * size0_y;
   unsigned long ind1_start_x = (offset_min_x) * size1_y;
   for (unsigned long n = 0; n < m._size; n++) {
      /* compute range of offset_y */
      unsigned long offset_min_y =
         ((pos_y + 1) > size0_y) ? (pos_y + 1 - size0_y) : 0;
      unsigned long offset_max_y =
         (pos_y < size1_y) ? pos_y : (size1_y - 1);
      unsigned long offset_range_y = offset_max_y - offset_min_y;
      /* initialize indices */
      unsigned long ind0 = ind0_start_x + (pos_y - offset_min_y);
      unsigned long ind1 = ind1_start_x + offset_min_y;
      for (unsigned long o_x = offset_min_x; o_x <= offset_max_x; o_x++) {
         for (unsigned long o_y = offset_min_y; o_y < offset_max_y; o_y++) {
            /* update result value */
            m._data[n] += m0._data[ind0] * m1._data[ind1];
            /* update linear positions */
            ind0--;
            ind1++;
         }
         /* update last result value */
         m._data[n] += m0._data[ind0] * m1._data[ind1];
         /* update linear positions */
         ind0 = ind0 + offset_range_y - size0_y;
         ind1 = ind1 - offset_range_y + size1_y;
      }
      /* update position */
      pos_y++;
      if (pos_y == pos_bound_y) {
         /* reset y position, increment x position */
         pos_y = pos_start_y;
         pos_x++;
         /* update range of offset_x */
         offset_min_x = ((pos_x + 1) > size0_x) ? (pos_x + 1 - size0_x) : 0;
         offset_max_x = (pos_x < size1_x) ? pos_x : (size1_x - 1);
         ind0_start_x = (pos_x - offset_min_x) * size0_y;
         ind1_start_x = (offset_min_x) * size1_y;
      }
   }
   return m;
}

/*
 * Compute convolution (for matrices of arbitrary dimensionality).
 */
template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::compute_conv(
   const matrix<T,Syn>& m0,
   const matrix<T,Syn>& m1,
   bool is_cropped,
   bool is_strict)
{
   /* get number of dimensions */
   unsigned long n_dims0 = m0._dims.size();
   unsigned long n_dims1 = m1._dims.size();
   unsigned long n_dims = (n_dims0 > n_dims1) ? n_dims0 : n_dims1;
   /* call faster versions for 1D matrices (if possible) */
   if (n_dims == 1) {
      return matrix<T,Syn>::compute_conv_1D(
         m0, m1, is_cropped, is_strict
      );
   } 
   /* extend dimension arrays */
   array<unsigned long> dims0 = matrix<T,Syn>::extend_dims(m0._dims, n_dims);
   array<unsigned long> dims1 = matrix<T,Syn>::extend_dims(m1._dims, n_dims);
   /* call faster versions for 2D matrices (if possible) */
   if (n_dims == 2) {
      return matrix<T,Syn>::compute_conv_2D(
         m0, m1, dims0, dims1, is_cropped, is_strict
      );
   }
   /* compute dimensions of resulting matrix and also     */ 
   /* compute range of position within full result matrix */
   array<unsigned long> dims(n_dims);
   array<unsigned long> pos_start(n_dims);
   array<unsigned long> pos_bound(n_dims);
   for (unsigned long d = 0; d < n_dims; d++) {
      if (!is_cropped) {
         /* set dimensions for full result matrix (start position is zero) */
         dims[d] = ((dims0[d] > 0) && (dims1[d] > 0)) ?
            (dims0[d] + dims1[d] - 1) : 0;
      } else if (!is_strict) {
         /* set dimensions for result matrix no larger than left input */
         dims[d] = ((dims0[d] > 0) && (dims1[d] > 0)) ? (dims0[d]) : 0;
         /* set start position for result matrix no larger than left input */
         pos_start[d] = dims1[d]/2;
      } else {
         /* set dimensions for central portion of result matrix */
         dims[d] = ((dims0[d] >= dims1[d]) && (dims1[d] > 0)) ?
            (dims0[d] - dims1[d] + 1) : 0;
         /* set start position for central portion of result matrix */
         pos_start[d] = (dims1[d] > 0) ? (dims1[d] - 1) : 0;
      }
      /* set bound on position range */
      pos_bound[d] = pos_start[d] + dims[d];
   }
   /* initialize result */
   matrix<T,Syn> m(dims);
   if (m._size == 0)
      return m;
   /* compute step sizes along each dimension */
   array<unsigned long> sizes0 = matrix<T,Syn>::dims_sizes(dims0);
   array<unsigned long> sizes1 = matrix<T,Syn>::dims_sizes(dims1);
   /* allocate storage for linear index range along each dimension */
   array<unsigned long> ind0_start(n_dims);
   array<unsigned long> ind1_start(n_dims);
   array<unsigned long> ind_cnt_max(n_dims);
   array<unsigned long> ind_cnt(n_dims);
   /* compute initial linear index range along each dimension */
   for (unsigned long d = 0; d < n_dims; d++) {
      unsigned long p = pos_start[d];
      unsigned long offset_min = ((p + 1) > dims0[d]) ? (p + 1 - dims0[d]) : 0;
      unsigned long offset_max = (p < dims1[d]) ? p : (dims1[d] - 1);
      ind0_start[d] = (p - offset_min) * sizes0[d];
      ind1_start[d] = (offset_min) * sizes1[d];
      ind_cnt_max[d] = offset_max - offset_min;
   }
   /* initialize position in result */
   array<unsigned long> pos(pos_start);
   /* compute convolution */   
   for (unsigned long n = 0; n < m._size; n++) {
      /* initialize linear positions in m0 and m1 */
      unsigned long ind0 = 0;
      unsigned long ind1 = 0;
      for (unsigned long d = 0; d < n_dims; d++) {
         ind0 += ind0_start[d];
         ind1 += ind1_start[d];
         ind_cnt[d] = 0;
      }
      /* multiply and add corresponding elements */
      while (true) {
         /* update result value */
         m._data[n] += m0._data[ind0] * m1._data[ind1];
         /* update offset */
         unsigned long d = n_dims - 1;
         while ((ind_cnt[d] == ind_cnt_max[d]) && (d > 0)) {
            ind0 += ind_cnt_max[d] * sizes0[d];
            ind1 -= ind_cnt_max[d] * sizes1[d];
            ind_cnt[d] = 0;
            d--;
         }
         if (ind_cnt[d] < ind_cnt_max[d]) {
            ind0 -= sizes0[d];
            ind1 += sizes1[d];
            ind_cnt[d]++;
         } else {
            break;
         }
      }
      /* update position */
      {
         unsigned long d = n_dims - 1;
         pos[d]++;
         while ((pos[d] == pos_bound[d]) && (d > 0)) {
            /* reset position along current dimension */
            unsigned long p = pos_start[d];
            pos[d] = p;
            /* reset linear index range along current dimension */
            unsigned long offset_min =
               ((p + 1) > dims0[d]) ? (p + 1 - dims0[d]) : 0;
            unsigned long offset_max =
               (p < dims1[d]) ? p : (dims1[d] - 1);
            ind0_start[d] = (p - offset_min) * sizes0[d];
            ind1_start[d] = (offset_min) * sizes1[d];
            ind_cnt_max[d] = offset_max - offset_min;      
            /* increment position along previous dimension */
            d--;
            pos[d]++;
         }
         /* update linear index range along current dimension */
         unsigned long p = pos[d];
         unsigned long offset_min =
            ((p + 1) > dims0[d]) ? (p + 1 - dims0[d]) : 0;
         unsigned long offset_max =
            (p < dims1[d]) ? p : (dims1[d] - 1);
         ind0_start[d] = (p - offset_min) * sizes0[d];
         ind1_start[d] = (offset_min) * sizes1[d];
         ind_cnt_max[d] = offset_max - offset_min; 
      }
   }
   return m;
}

/*
 * Dot product.
 */
template <typename T, typename Syn>
T dot(const matrix<T,Syn>& m0, const matrix<T,Syn>& m1) {
   auto_read_read_lock<const Syn> rrlock(m0, m1);
   matrix<T,Syn>::assert_dims_equal(m0._dims, m1._dims);
   T t = T();
   for (unsigned long n = 0; n < m0._size; n++)
      t += (m0._data[n] * m1._data[n]);
   return t;
}

/*
 * Element product.
 * Compute the product of corresponding matrix elements.
 */
template <typename T, typename Syn>
matrix<T,Syn> prod(
   const matrix<T,Syn>& m0, const matrix<T,Syn>& m1)
{
   auto_read_read_lock<const synchronizable> rrlock(m0, m1);
   matrix<T,Syn>::assert_dims_equal(m0._dims, m1._dims);
   matrix<T,Syn> m(m0._dims);
   for (unsigned long n = 0; n < m0._size; n++)
      m._data[n] = (m0._data[n] * m1._data[n]);
   return m;
}

/***************************************************************************
 * Minimum, maximum, sum, product.
 ***************************************************************************/

/*
 * Minimum, maximum, sum, product of all matrix elements.
 */
template <typename T, typename Syn>
T min(const matrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   /* check argument */
   if (m._size == 0)
      throw ex_invalid_argument(
         "cannot compute min of empty matrix"
      );
   /* compute min */
   T t = m._data[0];
   for (unsigned long n = 1; n < m._size; n++) {
      if (m._data[n] < t)
         t = m._data[n];
   }
   return t;
}

template <typename T, typename Syn>
T max(const matrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   /* check argument */
   if (m._size == 0)
      throw ex_invalid_argument(
         "cannot compute max of empty matrix"
      );
   /* compute max */
   T t = m._data[0];
   for (unsigned long n = 1; n < m._size; n++) {
      if (m._data[n] > t)
         t = m._data[n];
   }
   return t;
}

template <typename T, typename Syn>
T sum(const matrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   T t = T();
   for (unsigned long n = 0; n < m._size; n++)
      t += m._data[n];
   return t;
}

template <typename T, typename Syn>
T prod(const matrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   T t(1);
   for (unsigned long n = 0; n < m._size; n++)
      t *= m._data[n];
   return t;
}

/*
 * Minimum, maximum, sum, product along specified dimension.
 */
template <typename T, typename Syn>
matrix<T,Syn> min(const matrix<T,Syn>& m, unsigned long d) {
   auto_read_lock<const Syn> rlock(m);
   /* get # steps along dimension */
   unsigned long n_dims  = m._dims.size();
   unsigned long n_steps = (d < n_dims) ? m._dims[d] : 1;
   /* check if operating along zero or singleton dimension */
   if ((n_steps <= 1) || (m._size == 0))
      return matrix<T,Syn>::copy(m);
   /* compute indices along first slice of dimension */
   array<unsigned long> inds = matrix<T,Syn>::linear_indices_slice(m._dims, d);
   /* compute step size along dimension */
   unsigned long step_size = matrix<T,Syn>::dims_size(m._dims, d); 
   /* compute dimensions of result matrix */
   array<unsigned long> dims_result(m._dims);
   dims_result[d] = 1;
   /* compute result matrix */
   matrix<T,Syn> m_result(dims_result);
   for (unsigned long n = 0; n < m_result._size; n++) {
      unsigned long ind = inds[n];
      m_result._data[n] = m._data[ind];
      for (unsigned long n_step = 1; n_step < n_steps; n_step++) {
         ind += step_size;
         if (m._data[ind] < m_result._data[n])
            m_result._data[n] = m._data[ind];
      }
   }
   return m_result;
}

template <typename T, typename Syn>
matrix<T,Syn> max(const matrix<T,Syn>& m, unsigned long d) {
   auto_read_lock<const Syn> rlock(m);
   /* get # steps along dimension */
   unsigned long n_dims  = m._dims.size();
   unsigned long n_steps = (d < n_dims) ? m._dims[d] : 1;
   /* check if operating along zero or singleton dimension */
   if ((n_steps <= 1) || (m._size == 0))
      return matrix<T,Syn>::copy(m);
   /* compute indices along first slice of dimension */
   array<unsigned long> inds = matrix<T,Syn>::linear_indices_slice(m._dims, d);
   /* compute step size along dimension */
   unsigned long step_size = matrix<T,Syn>::dims_size(m._dims, d); 
   /* compute dimensions of result matrix */
   array<unsigned long> dims_result(m._dims);
   dims_result[d] = 1;
   /* compute result matrix */
   matrix<T,Syn> m_result(dims_result);
   for (unsigned long n = 0; n < m_result._size; n++) {
      unsigned long ind = inds[n];
      m_result._data[n] = m._data[ind];
      for (unsigned long n_step = 1; n_step < n_steps; n_step++) {
         ind += step_size;
         if (m._data[ind] > m_result._data[n])
            m_result._data[n] = m._data[ind];
      }
   }
   return m_result;
}

template <typename T, typename Syn>
matrix<T,Syn> sum(const matrix<T,Syn>& m, unsigned long d) {
   auto_read_lock<const Syn> rlock(m);
   /* get # steps along dimension */
   unsigned long n_dims  = m._dims.size();
   unsigned long n_steps = (d < n_dims) ? m._dims[d] : 1;
   /* check if operating along singleton dimension */
   if (n_steps == 1)
      return matrix<T,Syn>::copy(m);
   /* compute dimensions of result matrix */
   array<unsigned long> dims_result(m._dims);
   dims_result[d] = 1;
   /* initialize result matrix */
   matrix<T,Syn> m_result = matrix<T,Syn>::zeros(dims_result);
   /* check that result is nontrivial */
   if ((m_result._size > 0) && (m._size > 0)) {
      /* compute indices along first slice of dimension */
      array<unsigned long> inds = matrix<T,Syn>::linear_indices_slice(
         m._dims, d
      );
      /* compute step size along dimension */
      unsigned long step_size = matrix<T,Syn>::dims_size(m._dims, d);
      /* compute result matrix */
      for (unsigned long n = 0; n < m_result._size; n++) {
         unsigned long ind = inds[n];
         for (unsigned long n_step = 0; n_step < n_steps; n_step++) {
            m_result._data[n] += m._data[ind];
            ind += step_size;
         }
      }
   }
   return m_result;
}

template <typename T, typename Syn>
matrix<T,Syn> prod(const matrix<T,Syn>& m, unsigned long d) {
   auto_read_lock<const Syn> rlock(m);
   /* get # steps along dimension */
   unsigned long n_dims  = m._dims.size();
   unsigned long n_steps = (d < n_dims) ? m._dims[d] : 1;
   /* check if operating along singleton dimension */
   if (n_steps == 1)
      return matrix<T,Syn>::copy(m);
   /* compute dimensions of result matrix */
   array<unsigned long> dims_result(m._dims);
   dims_result[d] = 1;
   /* initialize result matrix */
   matrix<T,Syn> m_result = matrix<T,Syn>::ones(dims_result);
   /* check that result is nontrivial */
   if ((m_result._size > 0) && (m._size > 0)) {
      /* compute indices along first slice of dimension */
      array<unsigned long> inds = matrix<T,Syn>::linear_indices_slice(
         m._dims, d
      );
      /* compute step size along dimension */
      unsigned long step_size = matrix<T,Syn>::dims_size(m._dims, d);
      /* compute result matrix */
      for (unsigned long n = 0; n < m_result._size; n++) {
         unsigned long ind = inds[n];
         for (unsigned long n_step = 0; n_step < n_steps; n_step++) {
            m_result._data[n] *= m._data[ind];
            ind += step_size;
         }
      }
   }
   return m_result;
}

/***************************************************************************
 * Cumulative sum, product.
 ***************************************************************************/

/*
 * Cumulative sum, product of all matrix elements.
 */
template <typename T, typename Syn>
matrix<T,Syn> cumsum(const matrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   matrix<T,Syn> m_sum(m._dims);
   T t = T();
   for (unsigned long n = 0; n < m._size; n++) {
      t += m._data[n];
      m_sum._data[n] = t;
   }
   return m_sum;
}

template <typename T, typename Syn>
matrix<T,Syn> cumprod(const matrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   matrix<T,Syn> m_prod(m._dims);
   T t(1);
   for (unsigned long n = 0; n < m._size; n++) {
      t *= m._data[n];
      m_prod._data[n] = t;
   }
   return m_prod;
}

/*
 * Cumulative sum, product along specified dimension.
 */
template <typename T, typename Syn>
matrix<T,Syn> cumsum(const matrix<T,Syn>& m, unsigned long d) {
   auto_read_lock<const Syn> rlock(m);
   /* get # steps along dimension */
   unsigned long n_dims  = m._dims.size();
   unsigned long n_steps = (d < n_dims) ? m._dims[d] : 1;
   /* check if operating along singleton dimension */
   if (n_steps == 1)
      return matrix<T,Syn>::copy(m);
   /* initialize result matrix */
   matrix<T,Syn> m_result(m._dims);
   /* check that result is nontrivial */
   if (m_result._size > 0) {
      /* compute indices along first slice of dimension */
      array<unsigned long> inds = matrix<T,Syn>::linear_indices_slice(
         m._dims, d
      );
      /* compute step size along dimension */
      unsigned long step_size = matrix<T,Syn>::dims_size(m._dims, d);
      /* compute result matrix */
      unsigned long n_inds = inds.size();
      for (unsigned long n = 0; n < n_inds; n++) {
         unsigned long ind = inds[n];
         T t = T();
         for (unsigned long n_step = 0; n_step < n_steps; n_step++) {
            t += m._data[ind];
            m_result._data[ind] = t;
            ind += step_size;
         }
      }
   }
   return m_result;
}

template <typename T, typename Syn>
matrix<T,Syn> cumprod(const matrix<T,Syn>& m, unsigned long d) {
   auto_read_lock<const Syn> rlock(m);
   /* get # steps along dimension */
   unsigned long n_dims  = m._dims.size();
   unsigned long n_steps = (d < n_dims) ? m._dims[d] : 1;
   /* check if operating along singleton dimension */
   if (n_steps == 1)
      return matrix<T,Syn>::copy(m);
   /* initialize result matrix */
   matrix<T,Syn> m_result(m._dims);
   /* check that result is nontrivial */
   if (m_result._size > 0) {
       /* compute indices along first slice of dimension */
      array<unsigned long> inds = matrix<T,Syn>::linear_indices_slice(
         m._dims, d
      );
      /* compute step size along dimension */
      unsigned long step_size = matrix<T,Syn>::dims_size(m._dims, d);
      /* compute result matrix */
      unsigned long n_inds = inds.size();
      for (unsigned long n = 0; n < n_inds; n++) {
         unsigned long ind = inds[n];
         T t(1);
         for (unsigned long n_step = 0; n_step < n_steps; n_step++) {
            t *= m._data[ind];
            m_result._data[ind] = t;
            ind += step_size;
         }
      }
   }
   return m_result;
}

/***************************************************************************
 * Mean, variance.
 ***************************************************************************/

/*
 * Mean, variance of all matrix elements.
 */
template <typename T, typename Syn>
T mean(const matrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   /* check argument */
   if (m._size == 0)
      throw ex_invalid_argument(
         "cannot compute mean of empty matrix"
      );
   /* compute mean */
   T t_mean = T();
   for (unsigned long n = 0; n < m._size; n++)
      t_mean += m._data[n];
   t_mean /= T(m._size);
   return t_mean;
}

template <typename T, typename Syn>
T var(const matrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   /* check argument */
   if (m._size == 0)
      throw ex_invalid_argument(
         "cannot compute variance of empty matrix"
      );
   /* compute mean */
   T t_mean = T();
   for (unsigned long n = 0; n < m._size; n++)
      t_mean += m._data[n];
   t_mean /= T(m._size);
   /* compute variance */
   T t_var = T();
   for (unsigned long n = 0; n < m._size; n++) {
      T diff = m._data[n] - t_mean;
      t_var += diff * diff;
   }
   t_var /= T(m._size);
   return t_var;
}

/*
 * Mean, variance along specified dimension.
 */
template <typename T, typename Syn>
matrix<T,Syn> mean(const matrix<T,Syn>& m, unsigned long d) {
   auto_read_lock<const Syn> rlock(m);
   /* get # steps along dimension */
   unsigned long n_dims  = m._dims.size();
   unsigned long n_steps = (d < n_dims) ? m._dims[d] : 1;
   /* check if operating along singleton dimension */
   if (n_steps == 1)
      return matrix<T,Syn>::copy(m);
   /* compute dimensions of result matrix */
   array<unsigned long> dims_result(m._dims);
   dims_result[d] = 1;
   /* initialize result matrix */
   matrix<T,Syn> m_result = matrix<T,Syn>::zeros(dims_result);
   /* check that result is nontrivial */
   if ((m_result._size > 0) && (m._size > 0)) {
      /* compute indices along first slice of dimension */
      array<unsigned long> inds = matrix<T,Syn>::linear_indices_slice(
         m._dims, d
      );
      /* compute step size along dimension */
      unsigned long step_size = matrix<T,Syn>::dims_size(m._dims, d);
      /* compute result matrix */
      for (unsigned long n = 0; n < m_result._size; n++) {
         unsigned long ind = inds[n];
         for (unsigned long n_step = 0; n_step < n_steps; n_step++) {
            m_result._data[n] += m._data[ind];
            ind += step_size;
         }
         m_result._data[n] /= T(n_steps);
      }
   }
   return m_result;
}

template <typename T, typename Syn>
matrix<T,Syn> var(const matrix<T,Syn>& m, unsigned long d) {
   auto_read_lock<const Syn> rlock(m);
   /* get # steps along dimension */
   unsigned long n_dims  = m._dims.size();
   unsigned long n_steps = (d < n_dims) ? m._dims[d] : 1;
   /* check if operating along singleton dimension */
   if (n_steps == 1)
      return matrix<T,Syn>::zeros(m._dims);
   /* compute dimensions of result matrix */
   array<unsigned long> dims_result(m._dims);
   dims_result[d] = 1;
   /* initialize result matrix */
   matrix<T,Syn> m_result = matrix<T,Syn>::zeros(dims_result);
   /* check that result is nontrivial */
   if ((m_result._size > 0) && (m._size > 0)) {
      /* compute indices along first slice of dimension */
      array<unsigned long> inds = matrix<T,Syn>::linear_indices_slice(
         m._dims, d
      );
      /* compute step size along dimension */
      unsigned long step_size = matrix<T,Syn>::dims_size(m._dims, d);
      /* compute result matrix */
      for (unsigned long n = 0; n < m_result._size; n++) {
         /* compute mean */
         unsigned long ind = inds[n];
         T t_mean = T();
         for (unsigned long n_step = 0; n_step < n_steps; n_step++) {
            t_mean += m._data[ind];
            ind += step_size;
         }
         t_mean /= T(n_steps);
         /* compute variance */
         ind = inds[n];
         for (unsigned long n_step = 0; n_step < n_steps; n_step++) {
            T diff = m._data[ind] - t_mean;
            m_result._data[n] += diff * diff;
            ind += step_size;
         }
         m_result._data[n] /= T(n_steps);
      }
   }
   return m_result;
}

/***************************************************************************
 * Gradient.
 ***************************************************************************/

/*
 * Gradient along specified dimension.
 * Optionally specify the spacing (defaults to unit spacing) of elements
 * along the dimension.
 */
template <typename T, typename Syn>
matrix<T,Syn> gradient(const matrix<T,Syn>& m, unsigned long d) {
   return gradient(m, d, T(1));
}

template <typename T, typename Syn>
matrix<T,Syn> gradient(const matrix<T,Syn>& m, unsigned long d, const T& t) {
   const T t_space(t);
   auto_read_lock<const Syn> rlock(m);
   /* get # steps along dimension */
   unsigned long n_dims  = m._dims.size();
   unsigned long n_steps = (d < n_dims) ? m._dims[d] : 1;
   /* compute dimensions of gradient matrix */
   array<unsigned long> dims = matrix<T,Syn>::extend_dims(m._dims, d+1);
   dims[d] = (n_steps > 0) ? (n_steps-1) : 0;
   /* initialize gradient matrix */
   matrix<T,Syn> m_diff(dims);
   /* check that result is nontrivial */
   if (m_diff._size > 0) {
      /* compute indices along first slice of dimension */
      array<unsigned long> inds_m = matrix<T,Syn>::linear_indices_slice(
         m._dims, d
      );
      array<unsigned long> inds_m_diff = matrix<T,Syn>::linear_indices_slice(
         m_diff._dims, d
      );
      /* compute step size along dimension */
      unsigned long step_size = matrix<T,Syn>::dims_size(m._dims, d);
      /* compute gradient matrix */
      unsigned long n_inds = inds_m.size();
      for (unsigned long n = 0; n < n_inds; n++) {
         unsigned long ind_m_diff = inds_m_diff[n];
         unsigned long ind_m_prev = inds_m[n];
         for (unsigned long n_step = 1; n_step < n_steps; n_step++) {
            unsigned long ind_m_curr = ind_m_prev + step_size;
            m_diff._data[ind_m_diff] =
               (m._data[ind_m_curr] - m._data[ind_m_prev]) / (t_space);
            ind_m_prev = ind_m_curr;
            ind_m_diff += step_size;
         }
      }
   }
   return m_diff;
}

/***************************************************************************
 * Logical tests.
 ***************************************************************************/

/*
 * Logical AND of elements (true iff all element(s) are nonzero).
 */
template <typename T, typename Syn>
bool matrix<T,Syn>::all() const {
   auto_read_lock<const Syn> rlock(*this);
   const T t_zero = T();
   bool flag = true;
   for (unsigned long n = 0; ((n < this->_size) && (flag)); n++)
      flag = (this->_data[n] != t_zero);
   return flag;
}

/*
 * Logical OR of elements (true iff at least one element is nonzero).
 */
template <typename T, typename Syn>
bool matrix<T,Syn>::some() const {
   auto_read_lock<const Syn> rlock(*this);
   const T t_zero = T();
   bool flag = false;
   for (unsigned long n = 0; ((n < this->_size) && (!flag)); n++)
      flag = (this->_data[n] != t_zero);
   return flag;
}

/*
 * Logical NOR of elements (true iff all element(s) are zero).
 */
template <typename T, typename Syn>
bool matrix<T,Syn>::none() const {
   auto_read_lock<const Syn> rlock(*this);
   const T t_zero = T();
   bool flag = true;
   for (unsigned long n = 0; ((n < this->_size) && (flag)); n++)
      flag = (this->_data[n] == t_zero);
   return flag;
}
      
/*
 * Exactly one element true (true iff exactly one element is nonzero).
 */
template <typename T, typename Syn>
bool matrix<T,Syn>::one() const {
   auto_read_lock<const Syn> rlock(*this);
   const T t_zero = T();
   unsigned long nonzero_count = 0;
   for (unsigned long n = 0; ((n < this->_size) && (nonzero_count < 2)); n++) {
      if (this->_data[n] != t_zero)
         nonzero_count++;
   }
   return (nonzero_count == 1);
}

/***************************************************************************
 * Find.
 ***************************************************************************/

/*
 * Find locations of nonzero elements.
 */
template <typename T, typename Syn>
array<unsigned long> matrix<T,Syn>::find() const {
   auto_read_lock<const Syn> rlock(*this);
   const T t_zero = T();
   unsigned long nonzero_count = 0;
   for (unsigned long n = 0; n < this->_size; n++) {
      if (this->_data[n] != t_zero)
         nonzero_count++;
   }
   array<unsigned long> inds(nonzero_count);
   for (unsigned long n = 0, n_ind = 0; n_ind < nonzero_count; n++) {
      if (this->_data[n] != t_zero) {
         inds[n_ind] = n;
         n_ind++;
      }
   }
   return inds;
}

/*
 * Find locations of elements with the specified value.
 */
template <typename T, typename Syn>
array<unsigned long> matrix<T,Syn>::find_value(const T& t) const {
   auto_read_lock<const Syn> rlock(*this);
   const T t_copy(t);
   unsigned long count = 0;
   for (unsigned long n = 0; n < this->_size; n++) {
      if (this->_data[n] == t_copy)
         count++;
   }
   array<unsigned long> inds(count);
   for (unsigned long n = 0, n_ind = 0; n_ind < count; n++) {
      if (this->_data[n] == t_copy) {
         inds[n_ind] = n;
         n_ind++;
      }
   }
   return inds;
}

/***************************************************************************
 * Diagonal.
 ***************************************************************************/

/*
 * Diagonal entries (returns diagonal of matrix in a vector).
 */
template <typename T, typename Syn>
matrix<T,Syn> diag(const matrix<T,Syn>& m) {
   return m.diag();
}

template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::diag() const {
   auto_read_lock<const Syn> rlock(*this);
   /* compute increment size along diagonal */
   array<unsigned long> sizes = matrix<T,Syn>::dims_sizes(_dims);
   unsigned long n_dims = sizes.size();
   unsigned long step = 0;
   for (unsigned long n = 0; n < n_dims; n++)
      step += sizes[n];
   /* get length of diagonal */
   unsigned long n_dgnl = _dims[(matrix<T,Syn>::dims_min(_dims))];
   /* retreive diagonal elements */
   matrix<T,Syn> dgnl(n_dgnl, 1);
   for (unsigned long n = 0, pos = 0; n < n_dgnl; n++, pos += step)
      dgnl._data[n] = this->_data[pos];
   return dgnl;
}

/***************************************************************************
 * Inverse and row reduce.
 ***************************************************************************/

/*
 * Row reduce a 2D matrix (in place).
 *
 * Put the matrix in either row echelon form or reduced row echelon form.
 * Return the factor by which the determinant has been scaled.
 */
template <typename T, typename Syn>
T matrix<T,Syn>::row_reduce(const T& tol, bool to_rref) {
   /* check that matrix is 2D */
   unsigned long n_dims = _dims.size();
   if (n_dims != 2)
      throw ex_invalid_argument("matrix must be 2D");
   /* get matrix dimensions */
   unsigned long size_x = _dims[0];
   unsigned long size_y = _dims[1];
   unsigned long size_min = (size_x < size_y) ? size_x : size_y;
   /* define constants */
   const T t_zero = T();
   const T t_one(1);
   /* compute tolerance for matrix */
   T tolerance = tol;
   {
      /* compute infinity norm of matrix (maximum row sum) */
      T norm = T();
      for (unsigned long r = 0; r < this->_size; r += size_y) {
         T row_sum = T();
         for (unsigned long y = 0; y < size_y; y++)
            row_sum += abs(this->_data[r + y]);
         if (abs(row_sum) > abs(norm))
            norm = row_sum;
      }
      /* scale tolerance */
      tolerance *= T((size_x > size_y) ? size_x : size_y) * norm;
   }
   /* row reduce */
   T det_scale_factor = t_one;
   unsigned long x_ind = 0;
   for (unsigned long y = 0; y < size_min; y++) {
      /* find pivot in current column on or below current row */
      unsigned long pivot_x_ind = x_ind;
      T             pivot_val   = this->_data[pivot_x_ind + y];
      for (unsigned long r = x_ind + size_y; r < this->_size; r += size_y) {
         unsigned long ind = r + y;
         if (abs(this->_data[ind]) > abs(pivot_val)) {
            pivot_x_ind = r;
            pivot_val   = this->_data[ind];
         }
      }
      /* check that valid pivot was found */
      if (abs(pivot_val) > abs(tolerance)) {
         /* swap pivot row into current row */
         if (pivot_x_ind != x_ind) {
            for (unsigned long c = y; c < size_y; c++) {
               unsigned long ind       = x_ind + c;
               unsigned long pivot_ind = pivot_x_ind + c;
               T temp                 = this->_data[ind];
               this->_data[ind]       = this->_data[pivot_ind];
               this->_data[pivot_ind] = temp;
            }
            det_scale_factor = -det_scale_factor;
         }
         /* reduce portion of matrix above and on pivot row (if requested) */
         if (to_rref) { 
            /* divide pivot row by pivot value */
            T pivot_factor = t_one / pivot_val;
            this->_data[x_ind + y] = t_one;
            for (unsigned long c = y + 1; c < size_y; c++)
               this->_data[x_ind + c] *= pivot_factor;
            /* update pivot value and determinant scale factor */
            pivot_val = t_one;
            det_scale_factor *= pivot_factor;
            /* reduce rows above pivot row */
            for (unsigned long r = 0; r < x_ind; r += size_y) {
               T mult_factor = this->_data[r + y];
               this->_data[r + y] = t_zero;
               for (unsigned long c = y + 1; c < size_y; c++)
                  this->_data[r + c] -= mult_factor * this->_data[x_ind + c];
            }
         }
         /* reduce portion of matrix below pivot row */
         for (unsigned long r = x_ind + size_y; r < this->_size; r += size_y) {
            T mult_factor = this->_data[r + y] / pivot_val;
            this->_data[r + y] = t_zero;
            for (unsigned long c = y + 1; c < size_y; c++)
               this->_data[r + c] -= mult_factor * this->_data[x_ind + c];
         }
         /* increment current row */
         x_ind += size_y;
      } else {
         /* zero column on and below current row */
         for (unsigned long r = x_ind + y; r < this->_size; r += size_y)
            this->_data[r] = t_zero;
      }
   }
   return det_scale_factor;
}

/*
 * Inverse.
 */
template <typename T, typename Syn>
matrix<T,Syn> inv(const matrix<T,Syn>& m, const T& tol) {
   auto_read_lock<const Syn> rlock(m);
   return matrix<T,Syn>::compute_inv(m, tol);
}

template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::inv(const T& tol) {
   auto_write_lock<const Syn> wlock(*this);
   matrix<T,Syn> m = matrix<T,Syn>::compute_inv(*this, tol);
   matrix<T,Syn>::swap(*this, m);
   return *this;
}

template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::compute_inv(const matrix<T,Syn>& m, const T& tol) {
   /* check that matrix is 2D */
   if (m._dims.size() != 2)
      throw ex_invalid_argument("matrix must be 2D");
   /* check that matrix is square */
   if (m._dims[0] != m._dims[1])
      throw ex_invalid_argument("matrix must be square");
   /* define constants */
   const T t_zero = T();
   const T t_one(1);
   /* allocate augmented matrix */
   unsigned long N = m._dims[0];
   unsigned long two_N = N + N;
   matrix<T,Syn> m_eye(N, two_N);
   /* copy m */
   for (unsigned long dst = 0, src = 0; src < m._size; dst += N) {
      for (unsigned long y = 0; y < N; y++)
         m_eye._data[dst++] = m._data[src++];
   }
   /* copy identity */
   for (unsigned long dst = N; dst < m_eye._size; dst += (two_N+1))
      m_eye._data[dst] = t_one;
   /* row reduce */
   m_eye.row_reduce(tol, true);
   /* check that result is nonsingular */
   for (unsigned long n = 0; n < m_eye._size; n += (two_N+1)) {
      if (abs(m_eye._data[n]) == t_zero)
         throw ex_matrix_singular("cannot invert singular matrix");
   }
   /* grab inverse */
   matrix<T,Syn> m_inv(N, N);
   for (unsigned long src = N, dst = 0; dst < m_inv._size; src += N) {
      for (unsigned long y = 0; y < N; y++)
         m_inv._data[dst++] = m_eye._data[src++];
   }
   return m_inv;
}

/*
 * Row echelon form.
 */
template <typename T, typename Syn>
matrix<T,Syn> ref(const matrix<T,Syn>& m, const T& tol) {
   matrix<T,Syn> m_ref(m);
   m_ref.ref(tol);
   return m_ref;
}

template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::ref(const T& tol) {
   auto_write_lock<const Syn> wlock(*this);
   this->row_reduce(tol, false);
   return *this;
}

/*
 * Reduced row echelon form.
 */
template <typename T, typename Syn>
matrix<T,Syn> rref(const matrix<T,Syn>& m, const T& tol) {
   matrix<T,Syn> m_rref(m);
   m_rref.rref(tol);
   return m_rref;
}

template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::rref(const T& tol) {
   auto_write_lock<const Syn> wlock(*this);
   this->row_reduce(tol, true);
   return *this;
}

/***************************************************************************
 * Determinant and rank.
 ***************************************************************************/

/*
 * Determinant.
 */
template <typename T, typename Syn>
T matrix<T,Syn>::det(const T& tol) const {
   /* make working copy of matrix */
   matrix<T,Syn> m(*this);
   /* check that matrix is 2D */
   if (m._dims.size() != 2)
      throw ex_invalid_argument("matrix must be 2D");
   /* check that matrix is square */
   if (m._dims[0] != m._dims[1])
      throw ex_invalid_argument("matrix must be square");
   /* row reduce */
   T det_scale_factor = m.row_reduce(tol, false);
   /* compute determinant */
   T d = T(1) / det_scale_factor;
   unsigned long N = m._dims[0];
   for (unsigned long n = 0; n < m._size; n += (N+1))
      d *= m._data[n];
   return d;
}

/*
 * Rank.
 */
template <typename T, typename Syn>
unsigned long matrix<T,Syn>::rank(const T& tol) const {
   /* make working copy of matrix */
   matrix<T,Syn> m(*this);
   /* row reduce */
   m.row_reduce(tol, false);
   /* get matrix dimensions */
   unsigned long size_x = _dims[0];
   unsigned long size_y = _dims[1];
   unsigned long size_min = (size_x < size_y) ? size_x : size_y;
   /* compute rank */
   const T t_zero = T();
   unsigned long r = 0;
   for (unsigned long n = 0, y = 0; y < size_min; y++) {
      if (abs(m._data[n]) != t_zero) {
         n += size_y;
         r++;
      }
      n++;
   }
   return r;
}

/***************************************************************************
 * Dimensions and size.
 ***************************************************************************/

/*
 * Get dimensionality of matrix.
 */
template <typename T, typename Syn>
unsigned long matrix<T,Syn>::dimensionality() const {
   auto_read_lock<const Syn> rlock(*this);
   return _dims.size();
}

/*
 * Get dimensions of matrix.
 */
template <typename T, typename Syn>
array<unsigned long> matrix<T,Syn>::dimensions() const {
   auto_read_lock<const Syn> rlock(*this);
   return _dims;
}

/* 
 * Transpose.
 * Reverse dimension order.
 * The input matrix is regarded as being at least two-dimensional.
 */
template <typename T, typename Syn>
matrix<T,Syn> transpose(const matrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   return matrix<T,Syn>::compute_transpose(m);
}

template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::transpose() {
   auto_write_lock<const Syn> wlock(*this);
   matrix<T,Syn> m = matrix<T,Syn>::compute_transpose(*this);
   matrix<T,Syn>::swap(*this, m);
   return *this;
}

template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::compute_transpose(const matrix<T,Syn>& m) {
   /* create ordering that reverses dimension order */
   unsigned long n_dims = m._dims.size();
   n_dims = (n_dims < 2) ? 2 : n_dims;
   array<unsigned long> order(n_dims);
   for (unsigned long n = 0; n < n_dims; n++)
      order[n] = n_dims - n - 1;
   /* permute dimensions */
   return matrix<T,Syn>::compute_permute_dimensions(m, order);
}

/*
 * Dimension shifting.
 * Shift dimensions of the matrix to the left, wrapping the leading
 * dimensions around to the right.
 */
template <typename T, typename Syn>
matrix<T,Syn> shift_dimensions(
   const matrix<T,Syn>& m,
   unsigned long n_shift)
{
   auto_read_lock<const Syn> rlock(m);
   return matrix<T,Syn>::compute_shift_dimensions(m, n_shift);
}

template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::shift_dimensions(
   unsigned long n_shift)
{
   auto_write_lock<const Syn> wlock(*this);
   matrix<T,Syn> m = matrix<T,Syn>::compute_shift_dimensions(*this, n_shift);
   matrix<T,Syn>::swap(*this, m);
   return *this;
}

template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::compute_shift_dimensions(
   const matrix<T,Syn>& m,
   unsigned long n_shift)
{
   /* create ordering that shifts dimensions */
   unsigned long n_dims = m._dims.size();
   array<unsigned long> order(n_dims);
   n_shift = (n_dims > 0) ? (n_shift - n_dims*(n_shift/n_dims)) : 0;
   unsigned long n = 0;
   for (unsigned long d = n_shift; d < n_dims; d++, n++)
      order[n] = d;
   for (unsigned long d = 0; d < n_shift; d++, n++)
      order[n] = d;
   /* permute dimensions */
   return matrix<T,Syn>::compute_permute_dimensions(m, order);
}

/*
 * Dimension permutation.
 */
template <typename T, typename Syn>
matrix<T,Syn> permute_dimensions(
   const matrix<T,Syn>& m,
   const array<unsigned long>& order)
{
   auto_read_lock<const Syn> rlock(m);
   return matrix<T,Syn>::compute_permute_dimensions(m, order);
}

template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::permute_dimensions(
   const array<unsigned long>& order)
{
   auto_write_lock<const Syn> wlock(*this);
   matrix<T,Syn> m = matrix<T,Syn>::compute_permute_dimensions(*this, order);
   matrix<T,Syn>::swap(*this, m);
   return *this;
}

template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::compute_permute_dimensions(
   const matrix<T,Syn>& m,
   const array<unsigned long>& order)
{
   /* check that order is a permutation */
   unsigned long n_dims_order = order.size();
   array<bool> used(n_dims_order);
   for (unsigned long n = 0; n < n_dims_order; n++) {
      unsigned long d = order[n];
      if (d < n_dims_order) {
         if (!used[d]) {
            used[d] = true;
            continue;
         }
      }
      throw ex_invalid_argument(
         "invalid dimension permutation specified"
      );
   }
   /* extend dimensions array and compute sizes */
   unsigned long n_dims_m = m._dims.size();
   unsigned long n_dims = (n_dims_m > n_dims_order) ? n_dims_m : n_dims_order;
   array<unsigned long> dims  = matrix<T,Syn>::extend_dims(m._dims, n_dims);
   array<unsigned long> sizes = matrix<T,Syn>::dims_sizes(dims);
   /* permute dimensions and sizes */
   array<unsigned long> dims_permute(n_dims);
   array<unsigned long> sizes_permute(n_dims);
   for (unsigned long n = 0; n < n_dims_order; n++) {
      dims_permute[n]  = dims[order[n]];
      sizes_permute[n] = sizes[order[n]];
   }
   for (unsigned long n = n_dims_order; n < n_dims; n++) {
      dims_permute[n]  = dims[n];
      sizes_permute[n] = sizes[n];
   }
   /* initialize result matrix */
   matrix<T,Syn> m_perm(dims_permute);
   /* permute data order */
   if (m_perm._size > 0) {
      /* compute corresponding linear indices in original matrix */
      array<unsigned long> start(n_dims);
      array<unsigned long> end(n_dims);
      for (unsigned long n = 0; n < n_dims; n++)
         end[n] = dims_permute[n] - 1;
      array< array<unsigned long> > range = matrix<T,Syn>::range_indices(
         start, end
      );
      array<unsigned long> inds = matrix<T,Syn>::linear_indices_custom(
         sizes_permute, range
      );
      /* copy data */
      for (unsigned long n = 0; n < m_perm._size; n++)
         m_perm._data[n] = m._data[inds[n]];
   }
   return m_perm;
}

/*
 * Singleton dimension elimination.
 * Remove singleton dimensions, leaving the matrix elements unchanged.
 */
template <typename T, typename Syn>
matrix<T,Syn> squeeze_dimensions(const matrix<T,Syn>& m) {
   matrix<T,Syn> m_squeeze(m);
   m_squeeze.squeeze_dimensions();
   return m_squeeze;
}

template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::squeeze_dimensions() {
   auto_write_lock<const Syn> wlock(*this);
   unsigned long n_dims = _dims.size();
   if (n_dims > 1) {
      /* count number of singleton dimensions */
      unsigned long n_dims_single = 0;
      for (unsigned long d = 0; d < n_dims; d++) {
         if (_dims[d] == 1)
            n_dims_single++;
      }
      /* check if all dimensions are singleton */
      unsigned long n_dims_squeeze = n_dims - n_dims_single;
      if (n_dims_squeeze == 0) {
         /* leave one singleton dimension */
         array<unsigned long> dims_squeeze(1, 1);
         _dims = dims_squeeze;
      } else {
         /* remove all singleton dimensions */
         array<unsigned long> dims_squeeze(n_dims_squeeze);
         unsigned long d_squeeze = 0;
         for (unsigned long d = 0; d < n_dims; d++) {
            if (_dims[d] != 1) {
               dims_squeeze[d_squeeze] = _dims[d];
               d_squeeze++;
            }
         }
         _dims = dims_squeeze;
      }
   }
   return *this;
}

/*
 * Get size along specific dimension.
 */
template <typename T, typename Syn>
unsigned long matrix<T,Syn>::size(unsigned long d) const {
   auto_read_lock<const Syn> rlock(*this);
   return _dims[d];
}

/***************************************************************************
 * Resize.
 ***************************************************************************/

/*
 * Resize (vector).
 */
template <typename T, typename Syn>
matrix<T,Syn> resize(
   const matrix<T,Syn>& m, unsigned long length)
{
   auto_read_lock<const Syn> rlock(m);
   return matrix<T,Syn>::compute_resize(m, length);
}
   
template <typename T, typename Syn>
matrix<T,Syn> resize(
   const matrix<T,Syn>& m, unsigned long length, const T& t)
{
   auto_read_lock<const Syn> rlock(m);
   return matrix<T,Syn>::compute_resize(m, length, t);
}

template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::resize(unsigned long length) {
   auto_write_lock<const Syn> wlock(*this);
   matrix<T,Syn> m = matrix<T,Syn>::compute_resize(*this, length);
   matrix<T,Syn>::swap(*this, m);
   return *this;
}

template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::resize(unsigned long length, const T& t) {
   auto_write_lock<const Syn> wlock(*this);
   matrix<T,Syn> m = matrix<T,Syn>::compute_resize(*this, length, t);
   matrix<T,Syn>::swap(*this, m);
   return *this;
}

template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::compute_resize(
   const matrix<T,Syn>& m, unsigned long length)
{
   array<unsigned long> dims(1, length);
   return matrix<T,Syn>::compute_resize(m, dims);
}

template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::compute_resize(
   const matrix<T,Syn>& m, unsigned long length, const T& t)
{
   array<unsigned long> dims(1, length);
   return matrix<T,Syn>::compute_resize(m, dims, t);
}

/*
 * Resize (2D matrix).
 */
template <typename T, typename Syn>
matrix<T,Syn> resize(
   const matrix<T,Syn>& m, unsigned long M, unsigned long N)
{
   auto_read_lock<const Syn> rlock(m);
   return matrix<T,Syn>::compute_resize(m, M, N);
}

template <typename T, typename Syn>
matrix<T,Syn> resize(
   const matrix<T,Syn>& m, unsigned long M, unsigned long N, const T& t)
{
   auto_read_lock<const Syn> rlock(m);
   return matrix<T,Syn>::compute_resize(m, M, N, t);
}

template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::resize(
   unsigned long M, unsigned long N)
{
   auto_write_lock<const Syn> wlock(*this);
   matrix<T,Syn> m = matrix<T,Syn>::compute_resize(*this, M, N);
   matrix<T,Syn>::swap(*this, m);
   return *this;
}

template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::resize(
   unsigned long M, unsigned long N, const T& t)
{
   auto_write_lock<const Syn> wlock(*this);
   matrix<T,Syn> m = matrix<T,Syn>::compute_resize(*this, M, N, t);
   matrix<T,Syn>::swap(*this, m);
   return *this;
}

template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::compute_resize(
   const matrix<T,Syn>& m, unsigned long M, unsigned long N)
{
   array<unsigned long> dims(2);
   dims[0] = M;
   dims[1] = N;
   return matrix<T,Syn>::compute_resize(m, dims);
}

template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::compute_resize(
   const matrix<T,Syn>& m, unsigned long M, unsigned long N, const T& t)
{
   array<unsigned long> dims(2);
   dims[0] = M;
   dims[1] = N;
   return matrix<T,Syn>::compute_resize(m, dims, t);
}

/*
 * Resize (multi-dimensional matrix).
 * The dimensionality of the matrix must not change.
 */
template <typename T, typename Syn>
matrix<T,Syn> resize(
   const matrix<T,Syn>& m, const array<unsigned long>& dims)
{
   auto_read_lock<const Syn> rlock(m);
   return matrix<T,Syn>::compute_resize(m, dims);
}

template <typename T, typename Syn>
matrix<T,Syn> resize(
   const matrix<T,Syn>& m, const array<unsigned long>& dims, const T& t)
{
   auto_read_lock<const Syn> rlock(m);
   return matrix<T,Syn>::compute_resize(m, dims, t);
}

template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::resize(
   const array<unsigned long>& dims)
{
   auto_write_lock<const Syn> wlock(*this);
   matrix<T,Syn> m = matrix<T,Syn>::compute_resize(*this, dims);
   matrix<T,Syn>::swap(*this, m);
   return *this;
}

template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::resize(
   const array<unsigned long>& dims, const T& t)
{
   auto_write_lock<const Syn> wlock(*this);
   matrix<T,Syn> m = matrix<T,Syn>::compute_resize(*this, dims, t);
   matrix<T,Syn>::swap(*this, m);
   return *this;
}

template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::compute_resize(
   const matrix<T,Syn>& m, const array<unsigned long>& dims)
{
   auto_read_lock<const Syn> rlock(m);
   matrix<T,Syn> r(dims);
   matrix<T,Syn>::copy_overlap(m, r);
   return r;
}

template <typename T, typename Syn>
matrix<T,Syn> matrix<T,Syn>::compute_resize(
   const matrix<T,Syn>& m, const array<unsigned long>& dims, const T& t)
{
   auto_read_lock<const Syn> rlock(m);
   matrix<T,Syn> r(dims, t);
   matrix<T,Syn>::copy_overlap(m, r);
   return r;
}

/***************************************************************************
 * Reshape.
 ***************************************************************************/

/*
 * Reshape into a vector.
 */
template <typename T, typename Syn>
matrix<T,Syn> vector(const matrix<T,Syn>& m) {
   matrix<T,Syn> m_vector(m);
   m_vector.vector();
   return m_vector;
}

template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::vector() {
   auto_write_lock<const Syn> wlock(*this);
   _dims = array<unsigned long>(1, this->_size);
   return *this;
}

/*
 * Reshape into 2D M x N matrix.
 */
template <typename T, typename Syn>
matrix<T,Syn> reshape(
   const matrix<T,Syn>& m, unsigned long M, unsigned long N)
{
   matrix<T,Syn> m_reshaped(m);
   m_reshaped.reshape(M, N);
   return m_reshaped;
}

template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::reshape(unsigned long M, unsigned long N) {
   array<unsigned long> dims(2);
   dims[0] = M;
   dims[1] = N;
   return this->reshape(dims);
}

/*
 * Reshape into multi-dimensional matrix.
 */
template <typename T, typename Syn>
matrix<T,Syn> reshape(
   const matrix<T,Syn>& m, const array<unsigned long>& dims)
{
   matrix<T,Syn> m_reshaped(m);
   m_reshaped.reshape(dims);
   return m_reshaped;
}

template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::reshape(const array<unsigned long>& dims) {
   auto_write_lock<const Syn> wlock(*this);
   /* check argument */
   if (matrix<T,Syn>::dims_matrix_size(dims) != this->_size)
      throw ex_invalid_argument(
         "reshaping cannot change the size of a matrix"
      );
   /* reshape */
   _dims = dims;
   return *this;
}

/***************************************************************************
 * Replication.
 ***************************************************************************/

/*
 * Replicate matrix.
 */
template <typename T, typename Syn>
matrix<T,Syn> repmat(
   const matrix<T,Syn>& m, unsigned long M, unsigned long N)
{
   array<unsigned long> dims(2);
   dims[0] = M;
   dims[1] = N;
   return repmat(m, dims);
}
   
template <typename T, typename Syn>
matrix<T,Syn> repmat(
   const matrix<T,Syn>& m, const array<unsigned long>& dims)
{
   auto_read_lock<const Syn> rlock(m);
   /* compute resulting dimensions */
   unsigned long n_dims_m = m._dims.size();
   unsigned long n_dims_d = dims.size();
   unsigned long n_dims = (n_dims_m > n_dims_d) ? n_dims_m : n_dims_d;
   array<unsigned long> dims_m = matrix<T,Syn>::extend_dims(m._dims, n_dims);
   array<unsigned long> dims_d = matrix<T,Syn>::extend_dims(dims, n_dims);
   array<unsigned long> dims_rep(n_dims);
   for (unsigned long n = 0; n < n_dims; n++)
      dims_rep[n] = dims_m[n] * dims_d[n];
   /* initialize result */
   matrix<T,Syn> m_rep(dims_rep);
   /* replicate data */
   if (m_rep._size > 0) {
      /* compute indices of m within m_rep */
      array<unsigned long> start(n_dims);
      array<unsigned long> end(n_dims);
      for (unsigned long n = 0; n < n_dims; n++)
         end[n] = dims_m[n] - 1;
      array<unsigned long> inds_m = matrix<T,Syn>::linear_indices(
         dims_rep, start, end
      );
      /* compute offsets of replicated elements */
      array<unsigned long> step(n_dims);
      for (unsigned long n = 0; n < n_dims; n++) {
         step[n] = dims_m[n];
         end[n]  = dims_rep[n] - 1;
      }
      array<unsigned long> inds_offset = matrix<T,Syn>::linear_indices(
         dims_rep, start, step, end
      );
      /* copy data */
      unsigned long n_inds_offset = inds_offset.size();
      for (unsigned long n = 0; n < m._size; n++) {
         unsigned long ind = inds_m[n];
         for (unsigned long n_off = 0; n_off < n_inds_offset; n_off++) {
            unsigned long ind_off = inds_offset[n_off];
            m_rep._data[ind + ind_off] = m._data[n];
         }
      }
   }
   return m_rep;
}

/***************************************************************************
 * Concatenation.
 ***************************************************************************/

/*
 * Concatenate matrices.
 */
template <typename T, typename Syn>
matrix<T,Syn> vertcat(
   const matrix<T,Syn>& m_top, const matrix<T,Syn>& m_bottom)
{
   return concat(m_top, m_bottom, 0);
}

template <typename T, typename Syn>
matrix<T,Syn> horzcat(
   const matrix<T,Syn>& m_left, const matrix<T,Syn>& m_right)
{
   return concat(m_left, m_right, 1);
}

/*
 * Concatenate matrices along the specified dimension.
 */
template <typename T, typename Syn>
matrix<T,Syn> concat(
   const matrix<T,Syn>& m0, const matrix<T,Syn>& m1, unsigned long d)
{
   auto_read_read_lock<const Syn> rrlock(m0, m1);
   /* compute number of dimensions */
   unsigned long n_dims0 = m0._dims.size();
   unsigned long n_dims1 = m1._dims.size();
   unsigned long n_dims_max = (n_dims0 > n_dims1) ? n_dims0 : n_dims1;
   unsigned long n_dims = (d >= n_dims_max) ? (d + 1) : n_dims_max;
   /* extend dimension arrays */
   array<unsigned long> dims0 = matrix<T,Syn>::extend_dims(m0._dims, n_dims);
   array<unsigned long> dims1 = matrix<T,Syn>::extend_dims(m1._dims, n_dims);
   /* compute resulting dimensions */
   array<unsigned long> dims(n_dims);
   for (unsigned long n = 0; n < n_dims; n++) {
      if (n == d) {
         dims[n] = dims0[n] + dims1[n];
      } else if (dims0[n] != dims1[n]) {
         throw ex_matrix_dimension_mismatch(
            "matrix dimension mismatch during concat"
         );
      } else {
         dims[n] = dims0[n];
      }
   }
   /* initialize result */
   matrix<T,Syn> m(dims);
   /* copy m0 */
   if (m0._size > 0) {
      /* compute indices for placement of m0 */
      array<unsigned long> start(n_dims);
      array<unsigned long> end(n_dims);
      for (unsigned long n = 0; n < n_dims; n++)
         end[n] = dims0[n] - 1;
      array<unsigned long> inds = matrix<T,Syn>::linear_indices(
         dims, start, end
      );
      /* place m0 */
      for (unsigned long n = 0; n < m0._size; n++)
         m._data[inds[n]] = m0._data[n];
   }
   /* copy m0 */
   if (m1._size > 0) {
      /* compute indices for placement of m0 */
      array<unsigned long> start(n_dims);
      array<unsigned long> end(n_dims);
      for (unsigned long n = 0; n < n_dims; n++)
         end[n] = dims1[n] - 1;
      start[d] += dims0[d];
      end[d]   += dims0[d];
      array<unsigned long> inds = matrix<T,Syn>::linear_indices(
         dims, start, end
      );
      /* place m0 */
      for (unsigned long n = 0; n < m1._size; n++)
         m._data[inds[n]] = m1._data[n];
   }
   return m;
}

/***************************************************************************
 * Reverse.
 ***************************************************************************/

/*
 * Reverse element order (treating the matrix as a linear array).
 */
template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::reverse() {
   this->array<T,Syn>::reverse();
   return *this;
}

/*
 * Reverse element order along the specified dimension.
 */
template <typename T, typename Syn>
matrix<T,Syn>& matrix<T,Syn>::reverse(unsigned long d) {
   auto_write_lock<const Syn> wlock(*this);
   /* get # steps along dimension */
   unsigned long n_dims  = _dims.size();
   unsigned long n_steps = (d < n_dims) ? _dims[d] : 1;
   /* check if operating along zero or singleton dimension */
   if ((n_steps > 1) && (this->_size > 0)) {
      /* compute indices along first slice of dimension */
      array<unsigned long> inds = matrix<T,Syn>::linear_indices_slice(_dims, d);
      /* compute step size and max offset along dimension */
      unsigned long step_size  = matrix<T,Syn>::dims_size(_dims, d);
      unsigned long offset_max = (n_steps - 1) * step_size;
      /* reverse element order */
      unsigned long n_inds = inds.size();
      for (unsigned long n = 0; n < n_inds; n++) {
         unsigned long ind_low  = inds[n];
         unsigned long ind_high = ind_low + offset_max;
         while (ind_low < ind_high) {
            /* swap elements */
            T temp                = this->_data[ind_low];
            this->_data[ind_low]  = this->_data[ind_high];
            this->_data[ind_high] = temp;
            ind_low  += step_size;
            ind_high -= step_size;
         }
      }
   }
   return *this;
}

/***************************************************************************
 * Random permutation.
 ***************************************************************************/
 
/*
 * Randomly permute along the specified dimension.
 * Return an index matrix mapping resulting position --> original position.
 */
template <typename T, typename Syn>
matrix<unsigned long> matrix<T,Syn>::randperm(unsigned long d) {
   rand_gen_uniform<> r;
   return this->randperm(d, r);
}

template <typename T, typename Syn>
matrix<unsigned long> matrix<T,Syn>::randperm(unsigned long d, rand_gen<>& r) {
   auto_write_lock<const Syn> wlock(*this);
   /* get # steps along dimension */
   unsigned long n_dims  = _dims.size();
   unsigned long n_steps = (d < n_dims) ? _dims[d] : 1;
   /* initialize index matrix */
   matrix<unsigned long> m_idx(_dims);
   /* check if operating along zero or singleton dimension */
   if ((n_steps > 1) && (this->_size > 0)) {
      /* compute indices along first slice of dimension */
      array<unsigned long> inds = matrix<T,Syn>::linear_indices_slice(_dims, d);
      /* compute step size dimension */
      unsigned long step_size = matrix<T,Syn>::dims_size(_dims, d);
      /* perform permutation */
      unsigned long n_inds = inds.size();
      for (unsigned long n = 0; n < n_inds; n++) {
         /* grab elements */
         unsigned long ind = inds[n];
         matrix<T,Syn> m(n_steps, 1);
         for (unsigned long n_step = 0; n_step < n_steps; n_step++) {
            m._data[n_step] = this->_data[ind];
            ind += step_size;
         }
         /* permute elements */
         array<unsigned long> m_order = m.randperm(r);
         /* store permutation result */
         ind = inds[n];
         for (unsigned long n_step = 0; n_step < n_steps; n_step++) {
            this->_data[ind] = m._data[n_step];
            m_idx._data[ind] = m_order[n_step];
            ind += step_size;
         }
      }
   }
   return m_idx;
}

/***************************************************************************
 * Sorting.
 ***************************************************************************/

/*
 * Sort along specified dimension.
 * Elements are sorted in ascending order according to the given comparison
 * functor.
 */
template <typename T, typename Syn>
void matrix<T,Syn>::sort(
   unsigned long d, const comparable_functor<T>& f)
{
   auto_write_lock<const Syn> wlock(*this);
   /* get # steps along dimension */
   unsigned long n_dims  = _dims.size();
   unsigned long n_steps = (d < n_dims) ? _dims[d] : 1;
   /* check if operating along zero or singleton dimension */
   if ((n_steps > 1) && (this->_size > 0)) {
      /* compute indices along first slice of dimension */
      array<unsigned long> inds = matrix<T,Syn>::linear_indices_slice(_dims, d);
      /* compute step size dimension */
      unsigned long step_size = matrix<T,Syn>::dims_size(_dims, d);
      /* perform sort */
      unsigned long n_inds = inds.size();
      for (unsigned long n = 0; n < n_inds; n++) {
         /* grab elements */
         unsigned long ind = inds[n];
         matrix<T,Syn> m(n_steps, 1);
         for (unsigned long n_step = 0; n_step < n_steps; n_step++) {
            m._data[n_step] = this->_data[ind];
            ind += step_size;
         }
         /* sort elements */
         m.sort(f);
         /* store sort result */
         ind = inds[n];
         for (unsigned long n_step = 0; n_step < n_steps; n_step++) {
            this->_data[ind] = m._data[n_step];
            ind += step_size;
         }
      }
   }
}

/*
 * Sort along specified dimension.
 * Return an index matrix mapping sorted position --> original position.
 */
template <typename T, typename Syn>
matrix<unsigned long> matrix<T,Syn>::sort_idx(
   unsigned long d, const comparable_functor<T>& f)
{
   auto_write_lock<const Syn> wlock(*this);
   /* get # steps along dimension */
   unsigned long n_dims  = _dims.size();
   unsigned long n_steps = (d < n_dims) ? _dims[d] : 1;
   /* initialize index matrix */
   matrix<unsigned long> m_idx(_dims);
   /* check if operating along zero or singleton dimension */
   if ((n_steps > 1) && (this->_size > 0)) {
      /* compute indices along first slice of dimension */
      array<unsigned long> inds = matrix<T,Syn>::linear_indices_slice(_dims, d);
      /* compute step size dimension */
      unsigned long step_size = matrix<T,Syn>::dims_size(_dims, d);
      /* perform sort */
      unsigned long n_inds = inds.size();
      for (unsigned long n = 0; n < n_inds; n++) {
         /* grab elements */
         unsigned long ind = inds[n];
         matrix<T,Syn> m(n_steps, 1);
         for (unsigned long n_step = 0; n_step < n_steps; n_step++) {
            m._data[n_step] = this->_data[ind];
            ind += step_size;
         }
         /* sort elements */
         array<unsigned long> m_order = m.sort_idx(f);
         /* store sort result */
         ind = inds[n];
         for (unsigned long n_step = 0; n_step < n_steps; n_step++) {
            this->_data[ind] = m._data[n_step];
            m_idx._data[ind] = m_order[n_step];
            ind += step_size;
         }
      }
   }
   return m_idx;
}

/***************************************************************************
 * Uniqueness.
 ***************************************************************************/

/*
 * Remove duplicate elements (as defined by the given comparison functor)
 * and arrange the remaining unique elements in sorted order as a vector.
 */
template <typename T, typename Syn>
void matrix<T,Syn>::unique(const comparable_functor<T>& f) {
   auto_write_lock<const Syn> wlock(*this);
   this->make_unique(f);
   _dims = array<unsigned long>(1, this->_size);
}

/*
 * Remove duplicate elements (as defined by the given comparison functor)
 * and arrange the remaining unique elements in sorted order as a vector.
 * In addition, return an index array containing the original positions
 * of the unique elements.
 */
template <typename T, typename Syn>
array<unsigned long> matrix<T,Syn>::unique_idx(const comparable_functor<T>& f) {
   auto_write_lock<const Syn> wlock(*this);
   array<unsigned long> idx = this->make_unique_idx(f);
   _dims = array<unsigned long>(1, this->_size);
   return idx;
}

} /* namespace matrices */
} /* namespace math */

#endif
