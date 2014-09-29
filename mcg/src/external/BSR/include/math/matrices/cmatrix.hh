/*
 * Complex matrix (thread-safe).
 */
#ifndef MATH__MATRICES__CMATRIX_HH
#define MATH__MATRICES__CMATRIX_HH

#include "concurrent/threads/synchronization/locks/auto_read_lock.hh"
#include "concurrent/threads/synchronization/locks/auto_read_read_lock.hh"
#include "concurrent/threads/synchronization/locks/auto_write_lock.hh"
#include "concurrent/threads/synchronization/locks/auto_write_write_lock.hh"
#include "concurrent/threads/synchronization/synchronizables/unsynchronized.hh"
#include "io/serialization/serial_input_stream.hh"
#include "io/serialization/serial_output_stream.hh"
#include "io/serialization/serializers.hh"
#include "io/streams/ios.hh"
#include "io/streams/iomanip.hh"
#include "io/streams/ostream.hh"
#include "io/streams/ostringstream.hh"
#include "lang/array.hh"
#include "lang/exceptions/ex_invalid_argument.hh"
#include "math/complex.hh"
#include "math/math.hh"
#include "math/matrices/matrix.hh"
#include "math/random/generators/rand_gen.hh"
#include "math/random/generators/rand_gen_uniform.hh"

namespace math {
namespace matrices {
/*
 * Imports.
 */
using concurrent::threads::synchronization::locks::auto_read_lock;
using concurrent::threads::synchronization::locks::auto_read_read_lock;
using concurrent::threads::synchronization::locks::auto_write_lock;
using concurrent::threads::synchronization::locks::auto_write_write_lock;
using concurrent::threads::synchronization::synchronizables::unsynchronized;
using io::serialization::serial_input_stream;
using io::serialization::serial_output_stream;
using io::serialization::serializer;
using io::serialization::serializers;
using io::serialization::container_serializer;
using io::streams::ios;
using io::streams::ostream;
using io::streams::ostringstream;
using lang::array;
using lang::exceptions::ex_invalid_argument;
using math::complex;
using math::random::generators::rand_gen;
using math::random::generators::rand_gen_uniform;

/*
 * Import standard math functions.
 */
using math::acos;
using math::asin;
using math::acosh;
using math::atanh;
using math::sqrt;
using math::log;
using math::pow;
using math::atan2;
using math::abs;
using math::arg;
using math::conj;
using math::eps;

/*
 * Declare prototypes for cmatrix template friend functions.
 */
template <typename T, typename Syn>
ostream& operator<< (ostream&, const cmatrix<T,Syn>&);

template <typename T, typename Syn>
cmatrix<T,Syn> operator+(const cmatrix<T,Syn>&, const T&);
template <typename T, typename Syn>
cmatrix<T,Syn> operator-(const cmatrix<T,Syn>&, const T&);
template <typename T, typename Syn>
cmatrix<T,Syn> operator*(const cmatrix<T,Syn>&, const T&);
template <typename T, typename Syn>
cmatrix<T,Syn> operator/(const cmatrix<T,Syn>&, const T&);

template <typename T, typename Syn>
cmatrix<T,Syn> operator+(const T&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> operator-(const T&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> operator*(const T&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> operator/(const T&, const cmatrix<T,Syn>&);

template <typename T, typename Syn>
cmatrix<T,Syn> operator+(const matrix<T,Syn>&, const complex<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> operator-(const matrix<T,Syn>&, const complex<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> operator*(const matrix<T,Syn>&, const complex<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> operator/(const matrix<T,Syn>&, const complex<T,Syn>&);

template <typename T, typename Syn>
cmatrix<T,Syn> operator+(const complex<T,Syn>&, const matrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> operator-(const complex<T,Syn>&, const matrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> operator*(const complex<T,Syn>&, const matrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> operator/(const complex<T,Syn>&, const matrix<T,Syn>&);
   
template <typename T, typename Syn>
cmatrix<T,Syn> operator+(const cmatrix<T,Syn>&, const complex<T>&);
template <typename T, typename Syn>
cmatrix<T,Syn> operator-(const cmatrix<T,Syn>&, const complex<T>&);
template <typename T, typename Syn>
cmatrix<T,Syn> operator*(const cmatrix<T,Syn>&, const complex<T>&);
template <typename T, typename Syn>
cmatrix<T,Syn> operator/(const cmatrix<T,Syn>&, const complex<T>&);

template <typename T, typename Syn>
cmatrix<T,Syn> operator+(const complex<T>&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> operator-(const complex<T>&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> operator*(const complex<T>&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> operator/(const complex<T>&, const cmatrix<T,Syn>&);

template <typename T, typename Syn>
cmatrix<T,Syn> operator+(const cmatrix<T,Syn>&, const matrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> operator-(const cmatrix<T,Syn>&, const matrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> operator*(const cmatrix<T,Syn>&, const matrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> operator/(const cmatrix<T,Syn>&, const matrix<T,Syn>&);

template <typename T, typename Syn>
cmatrix<T,Syn> operator+(const matrix<T,Syn>&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> operator-(const matrix<T,Syn>&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> operator*(const matrix<T,Syn>&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> operator/(const matrix<T,Syn>&, const cmatrix<T,Syn>&);

template <typename T, typename Syn>
cmatrix<T,Syn> operator+(const cmatrix<T,Syn>&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> operator-(const cmatrix<T,Syn>&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> operator*(const cmatrix<T,Syn>&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> operator/(const cmatrix<T,Syn>&, const cmatrix<T,Syn>&);

template <typename T, typename Syn>
matrix<bool,Syn> operator==(const cmatrix<T,Syn>&, const T&);
template <typename T, typename Syn>
matrix<bool,Syn> operator!=(const cmatrix<T,Syn>&, const T&);
template <typename T, typename Syn>
matrix<bool,Syn> operator< (const cmatrix<T,Syn>&, const T&);
template <typename T, typename Syn>
matrix<bool,Syn> operator> (const cmatrix<T,Syn>&, const T&);
template <typename T, typename Syn>
matrix<bool,Syn> operator<=(const cmatrix<T,Syn>&, const T&);
template <typename T, typename Syn>
matrix<bool,Syn> operator>=(const cmatrix<T,Syn>&, const T&);

template <typename T, typename Syn>
matrix<bool,Syn> operator==(const T&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator!=(const T&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator< (const T&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator> (const T&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator<=(const T&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator>=(const T&, const cmatrix<T,Syn>&);

template <typename T, typename Syn>
matrix<bool,Syn> operator==(const cmatrix<T,Syn>&, const complex<T>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator!=(const cmatrix<T,Syn>&, const complex<T>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator< (const cmatrix<T,Syn>&, const complex<T>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator> (const cmatrix<T,Syn>&, const complex<T>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator<=(const cmatrix<T,Syn>&, const complex<T>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator>=(const cmatrix<T,Syn>&, const complex<T>&);

template <typename T, typename Syn>
matrix<bool,Syn> operator==(const complex<T>&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator!=(const complex<T>&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator< (const complex<T>&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator> (const complex<T>&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator<=(const complex<T>&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator>=(const complex<T>&, const cmatrix<T,Syn>&);

template <typename T, typename Syn>
matrix<bool,Syn> operator==(const cmatrix<T,Syn>&, const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator!=(const cmatrix<T,Syn>&, const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator< (const cmatrix<T,Syn>&, const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator> (const cmatrix<T,Syn>&, const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator<=(const cmatrix<T,Syn>&, const matrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator>=(const cmatrix<T,Syn>&, const matrix<T,Syn>&);

template <typename T, typename Syn>
matrix<bool,Syn> operator==(const matrix<T,Syn>&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator!=(const matrix<T,Syn>&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator< (const matrix<T,Syn>&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator> (const matrix<T,Syn>&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator<=(const matrix<T,Syn>&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator>=(const matrix<T,Syn>&, const cmatrix<T,Syn>&);

template <typename T, typename Syn>
matrix<bool,Syn> operator==(const cmatrix<T,Syn>&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator!=(const cmatrix<T,Syn>&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator< (const cmatrix<T,Syn>&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator> (const cmatrix<T,Syn>&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator<=(const cmatrix<T,Syn>&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
matrix<bool,Syn> operator>=(const cmatrix<T,Syn>&, const cmatrix<T,Syn>&);

template <typename T, typename Syn>
cmatrix<T,Syn> cos(const cmatrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> sin(const cmatrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> tan(const cmatrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> cosh(const cmatrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> sinh(const cmatrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> tanh(const cmatrix<T,Syn>&);

template <typename T, typename Syn>
cmatrix<T,Syn> acos(const cmatrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> asin(const cmatrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> atan(const cmatrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> acosh(const cmatrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> asinh(const cmatrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> atanh(const cmatrix<T,Syn>&);

template <typename T, typename Syn>
cmatrix<T,Syn> sqrt(const cmatrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> log(const cmatrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> exp(const cmatrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> pow(const cmatrix<T,Syn>&, int);
template <typename T, typename Syn>
cmatrix<T,Syn> pow(const cmatrix<T,Syn>&, const T&);
template <typename T, typename Syn>
cmatrix<T,Syn> pow(const cmatrix<T,Syn>&, const complex<T>&);

template <typename T, typename Syn>
matrix<T,Syn> atan2(const cmatrix<T,Syn>&, const cmatrix<T,Syn>&);

template <typename T, typename Syn>
matrix<T,Syn> abs(const cmatrix<T,Syn>&);
template <typename T, typename Syn>
matrix<T,Syn> arg(const cmatrix<T,Syn>&);

template <typename T, typename Syn>
cmatrix<T,Syn> conj(const cmatrix<T,Syn>&);

template <typename T, typename Syn>
cmatrix<T,Syn> conv(const cmatrix<T,Syn>&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> conv_crop(const cmatrix<T,Syn>&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> conv_crop_strict(const cmatrix<T,Syn>&, const cmatrix<T,Syn>&);

template <typename T, typename Syn>
complex<T> dot(const cmatrix<T,Syn>&, const cmatrix<T,Syn>&);

template <typename T, typename Syn>
cmatrix<T,Syn> prod(const cmatrix<T,Syn>&, const cmatrix<T,Syn>&);

template <typename T, typename Syn> complex<T> min(const cmatrix<T,Syn>&);
template <typename T, typename Syn> complex<T> max(const cmatrix<T,Syn>&);
template <typename T, typename Syn> complex<T> sum(const cmatrix<T,Syn>&);
template <typename T, typename Syn> complex<T> prod(const cmatrix<T,Syn>&);

template <typename T, typename Syn>
cmatrix<T,Syn> min(const cmatrix<T,Syn>&, unsigned long);
template <typename T, typename Syn>
cmatrix<T,Syn> max(const cmatrix<T,Syn>&, unsigned long);
template <typename T, typename Syn>
cmatrix<T,Syn> sum(const cmatrix<T,Syn>&, unsigned long);
template <typename T, typename Syn>
cmatrix<T,Syn> prod(const cmatrix<T,Syn>&, unsigned long);

template <typename T, typename Syn>
cmatrix<T,Syn> cumsum(const cmatrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> cumprod(const cmatrix<T,Syn>&);

template <typename T, typename Syn>
cmatrix<T,Syn> cumsum(const cmatrix<T,Syn>&, unsigned long);
template <typename T, typename Syn>
cmatrix<T,Syn> cumprod(const cmatrix<T,Syn>&, unsigned long);

template <typename T, typename Syn> complex<T> mean(const cmatrix<T,Syn>&);
template <typename T, typename Syn> complex<T> var(const cmatrix<T,Syn>&);

template <typename T, typename Syn>
cmatrix<T,Syn> mean(const cmatrix<T,Syn>&, unsigned long);
template <typename T, typename Syn>
cmatrix<T,Syn> var(const cmatrix<T,Syn>&, unsigned long);

template <typename T, typename Syn>
cmatrix<T,Syn> gradient(
   const cmatrix<T,Syn>&, unsigned long);
template <typename T, typename Syn>
cmatrix<T,Syn> gradient(
   const cmatrix<T,Syn>&, unsigned long, const complex<T>&);

template <typename T, typename Syn>
cmatrix<T,Syn> diag(const cmatrix<T,Syn>&);

template <typename T, typename Syn>
cmatrix<T,Syn> inv(const cmatrix<T,Syn>&, const T& = eps<T>());
template <typename T, typename Syn>
cmatrix<T,Syn> ref(const cmatrix<T,Syn>&, const T& = eps<T>());
template <typename T, typename Syn>
cmatrix<T,Syn> rref(const cmatrix<T,Syn>&, const T& = eps<T>());

template <typename T, typename Syn>
cmatrix<T,Syn> transpose(const cmatrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> shift_dimensions(
   const cmatrix<T,Syn>&, unsigned long);
template <typename T, typename Syn>
cmatrix<T,Syn> permute_dimensions(
   const cmatrix<T,Syn>&, const array<unsigned long>&);
template <typename T, typename Syn>
cmatrix<T,Syn> squeeze_dimensions(const cmatrix<T,Syn>&);

template <typename T, typename Syn>
cmatrix<T,Syn> resize(const cmatrix<T,Syn>&, unsigned long);
template <typename T, typename Syn>
cmatrix<T,Syn> resize(const cmatrix<T,Syn>&, unsigned long, const T&);
template <typename T, typename Syn>
cmatrix<T,Syn> resize(const cmatrix<T,Syn>&, unsigned long, const complex<T>&);

template <typename T, typename Syn>
cmatrix<T,Syn> resize(
   const cmatrix<T,Syn>&, unsigned long, unsigned long);
template <typename T, typename Syn>
cmatrix<T,Syn> resize(
   const cmatrix<T,Syn>&, unsigned long, unsigned long, const T&);
template <typename T, typename Syn>
cmatrix<T,Syn> resize(
   const cmatrix<T,Syn>&, unsigned long, unsigned long, const complex<T>&);

template <typename T, typename Syn>
cmatrix<T,Syn> resize(
   const cmatrix<T,Syn>&, const array<unsigned long>&);
template <typename T, typename Syn>
cmatrix<T,Syn> resize(
   const cmatrix<T,Syn>&, const array<unsigned long>&, const T&);
template <typename T, typename Syn>
cmatrix<T,Syn> resize(
   const cmatrix<T,Syn>&, const array<unsigned long>&, const complex<T>&);

template <typename T, typename Syn>
cmatrix<T,Syn> vector(const cmatrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> reshape(const cmatrix<T,Syn>&, unsigned long, unsigned long);
template <typename T, typename Syn>
cmatrix<T,Syn> reshape(const cmatrix<T,Syn>&, const array<unsigned long>&);

template <typename T, typename Syn>
cmatrix<T,Syn> repmat(const cmatrix<T,Syn>&, unsigned long, unsigned long);
template <typename T, typename Syn>
cmatrix<T,Syn> repmat(const cmatrix<T,Syn>&, const array<unsigned long>&);

template <typename T, typename Syn>
cmatrix<T,Syn> vertcat(const cmatrix<T,Syn>&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> horzcat(const cmatrix<T,Syn>&, const cmatrix<T,Syn>&);
template <typename T, typename Syn>
cmatrix<T,Syn> concat(
   const cmatrix<T,Syn>&, const cmatrix<T,Syn>&, unsigned long);

/*
 * Complex matrix class.
 */
template <typename T = double, typename Syn = unsynchronized>
class cmatrix : public matrix<complex<T>,Syn> {
public:
   /*
    * Friend classes.
    */
   template <typename U, typename S> friend class cmatrix;

   /*
    * Friend libraries.
    */
   friend class math::libraries::lib_image;
   friend class math::libraries::lib_matrix;
   friend class math::libraries::lib_signal;

   /*
    * Constructor.
    * Return the empty (0 x 0) matrix.
    */
   cmatrix();
   
   /*
    * Constructors for single element (1 x 1) matrix.
    */
   explicit cmatrix(const T&);
   explicit cmatrix(const complex<T>&);
   
   /*
    * Constructors for two-dimensional M x N matrices.
    */
   explicit cmatrix(
      unsigned long,                /* M */
      unsigned long                 /* N */
   );
 
   explicit cmatrix(
      unsigned long,                /* M */
      unsigned long,                /* N */
      const T&                      /* element initialization value */
   );
  
   explicit cmatrix(
      unsigned long,                /* M */
      unsigned long,                /* N */
      const complex<T>&             /* element initialization value */
   );
   
   /*
    * Constructors for multi-dimensional matrices.
    */
   explicit cmatrix(
      const array<unsigned long>&   /* dimensions */
   );
   
   explicit cmatrix(
      const array<unsigned long>&,  /* dimensions */
      const T&                      /* element initialization value */
   );
     
   explicit cmatrix(
      const array<unsigned long>&,  /* dimensions */
      const complex<T>&             /* element initialization value */
   );

   /*
    * Constructor from real-valued matrix.
    */
   explicit cmatrix(const matrix<T,Syn>&);
   
   /*
    * Constructor from real and imaginary components.
    */
   explicit cmatrix(
      const matrix<T,Syn>& /* real */, const matrix<T,Syn>& /* imag */
   );
   
   /*
    * Copy constructors.
    */
   cmatrix(const cmatrix<T,Syn>&);
   template <typename U, typename S> explicit cmatrix(const cmatrix<U,S>&);

   /*
    * Destructor.
    */
   virtual ~cmatrix();
   
   /*
    * Empty matrix (0 x 0).
    */
   static cmatrix<T,Syn> empty();
   
   /*
    * Identity matrix.
    * Create and return an N x N identity matrix.
    */
   static cmatrix<T,Syn> eye(unsigned long /* N */);
   
   /*
    * Zeros matrix.
    * Create and return an M x N matrix of zeros.
    */
   static cmatrix<T,Syn> zeros(
      unsigned long,                /* M */
      unsigned long                 /* N */
   );

   /*
    * Zeros matrix.
    * Create and return a matrix of zeros with the given dimensions.
    */
   static cmatrix<T,Syn> zeros(
      const array<unsigned long>&   /* dimensions */
   );
   
   /*
    * Ones matrix.
    * Create and return an M x N matrix of ones.
    */
   static cmatrix<T,Syn> ones(
      unsigned long,                /* M */
      unsigned long                 /* N */
   );
   
   /*
    * Ones matrix.
    * Create and return a matrix of ones with the given dimensions.
    */
   static cmatrix<T,Syn> ones(
      const array<unsigned long>&   /* dimensions */
   );

   /*
    * Random matrix.
    * Create and return an M x N matrix of random numbers drawn from the
    * given random number generator.  If no generator is specified, then
    * the matrix elements are uniformly distributed random numbers in the 
    * interval [0, 1] (inclusive).
    */ 
   static cmatrix<T,Syn> random(
      unsigned long,                /* M */
      unsigned long                 /* N */
   );

   static cmatrix<T,Syn> random(
      unsigned long,                /* M */
      unsigned long,                /* N */
      rand_gen<T>&                  /* generator */
   );

   static cmatrix<T,Syn> random(
      unsigned long,                /* M */
      unsigned long,                /* N */
      rand_gen< complex<T> >&       /* generator */
   );
  
   /*
    * Random matrix.
    * Create and return a matrix with the given dimensions in which each 
    * element is a random number drawn from the given random number generator.
    * If no generator is specified, then the matrix elements are uniformly 
    * distributed random numbers in the interval [0, 1] (inclusive).
    */
   static cmatrix<T,Syn> random(
      const array<unsigned long>&   /* dimensions */
   );

   static cmatrix<T,Syn> random(
      const array<unsigned long>&,  /* dimensions */
      rand_gen<T>&                  /* generator */
   );
   
   static cmatrix<T,Syn> random(
      const array<unsigned long>&,  /* dimensions */
      rand_gen< complex<T> >&       /* generator */
   );

   /*
    * Ramp.
    * Create and return a vector containing values from start to end at an 
    * increment of step.  The default step is one.
    */
   static cmatrix<T,Syn> ramp(
      const T&,                     /* start */
      const T&                      /* end */
   );
   
   static cmatrix<T,Syn> ramp(
      const T&,                     /* start */
      const T&,                     /* step */
      const T&                      /* end */
   );
   
   /*
    * Diagonal matrix.
    * Create and return a d-dimensional square matrix given its diagonal.
    */
   static cmatrix<T,Syn> diagonal(
      const matrix<T,Syn>&,         /* diagonal entries */
      unsigned long = 2             /* d (dimensionality) */
   );
      
   static cmatrix<T,Syn> diagonal(
      const cmatrix<T,Syn>&,        /* diagonal entries */
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
   static auto_ptr< cmatrix<T,Syn> > deserialize(
      serial_input_stream&,
      const serializer<T>& = serializers<T>::s_default()
   );

   /*
    * Element reference.
    */
   using matrix<complex<T>,Syn>::operator();

   /*
    * Direct element access.
    */
   using matrix<complex<T>,Syn>::data;

   /*
    * Subarray.
    * Return vector of elements in the specified linear index range.
    */
   cmatrix<T,Syn> subarray(
      unsigned long,                /* start index */ 
      unsigned long                 /* end index   */
   ) const;
   
   cmatrix<T,Syn> subarray(
      unsigned long,                /* start index */ 
      unsigned long,                /* step size   */
      unsigned long                 /* end index   */
   ) const;
   
   cmatrix<T,Syn> subarray(
      const array<unsigned long>&   /* indices of desired elements */
   ) const;
   
   /*
    * Submatrix.
    * Return matrix of elements in the specified range.
    */
   cmatrix<T,Syn> submatrix(
      const array<unsigned long>&,  /* start index along each dimension */
      const array<unsigned long>&   /* end index along each dimension */
   ) const;
   
   cmatrix<T,Syn> submatrix(
      const array<unsigned long>&,  /* start index along each dimension */
      const array<unsigned long>&,  /* step along each dimension */
      const array<unsigned long>&   /* end index along each dimension */
   ) const;
   
   cmatrix<T,Syn> submatrix(
      const array< array<unsigned long> >&   /* indices along each dimension */
   ) const;

   /*
    * Subarray assignment.
    * Assign value(s) to elements in the specified linear index range.
    */
   cmatrix<T,Syn>& subassign(
      unsigned long,                /* start index */ 
      unsigned long,                /* end index   */
      const T&                      /* value */
   );
   
   cmatrix<T,Syn>& subassign(
      unsigned long,                /* start index */ 
      unsigned long,                /* end index   */
      const complex<T>&             /* value */
   );
   
   cmatrix<T,Syn>& subassign(
      unsigned long,                /* start index */ 
      unsigned long,                /* end index   */
      const matrix<T,Syn>&          /* values */
   );
   
   cmatrix<T,Syn>& subassign(
      unsigned long,                /* start index */ 
      unsigned long,                /* end index   */
      const cmatrix<T,Syn>&         /* values */
   );

   cmatrix<T,Syn>& subassign(
      unsigned long,                /* start index */ 
      unsigned long,                /* step size   */
      unsigned long,                /* end index   */
      const T&                      /* value */
   );

   cmatrix<T,Syn>& subassign(
      unsigned long,                /* start index */ 
      unsigned long,                /* step size   */
      unsigned long,                /* end index   */
      const complex<T>&             /* value */
   );
   
   cmatrix<T,Syn>& subassign(
      unsigned long,                /* start index */ 
      unsigned long,                /* step size   */
      unsigned long,                /* end index   */
      const matrix<T,Syn>&          /* values */
   );

   cmatrix<T,Syn>& subassign(
      unsigned long,                /* start index */ 
      unsigned long,                /* step size   */
      unsigned long,                /* end index   */
      const cmatrix<T,Syn>&         /* values */
   );

   cmatrix<T,Syn>& subassign(
      const array<unsigned long>&,  /* indices of desired elements */
      const T&                      /* value */
   );
   
   cmatrix<T,Syn>& subassign(
      const array<unsigned long>&,  /* indices of desired elements */
      const complex<T>&             /* value */
   );
   
   cmatrix<T,Syn>& subassign(
      const array<unsigned long>&,  /* indices of desired elements */
      const matrix<T,Syn>&          /* values */
   );
   
   cmatrix<T,Syn>& subassign(
      const array<unsigned long>&,  /* indices of desired elements */
      const cmatrix<T,Syn>&         /* values */
   );

   /*
    * Submatrix assignment.
    * Assign value(s) to elements in the specified range.
    */
   cmatrix<T,Syn>& subassign(
      const array<unsigned long>&,  /* start index along each dimension */
      const array<unsigned long>&,  /* end index along each dimension */
      const T&                      /* value */
   );

   cmatrix<T,Syn>& subassign(
      const array<unsigned long>&,  /* start index along each dimension */
      const array<unsigned long>&,  /* end index along each dimension */
      const complex<T>&             /* value */
   );
   
   cmatrix<T,Syn>& subassign(
      const array<unsigned long>&,  /* start index along each dimension */
      const array<unsigned long>&,  /* end index along each dimension */
      const matrix<T,Syn>&          /* values */
   );
 
   cmatrix<T,Syn>& subassign(
      const array<unsigned long>&,  /* start index along each dimension */
      const array<unsigned long>&,  /* end index along each dimension */
      const cmatrix<T,Syn>&         /* values */
   );
   
   cmatrix<T,Syn>& subassign(
      const array<unsigned long>&,  /* start index along each dimension */
      const array<unsigned long>&,  /* step along each dimension */
      const array<unsigned long>&,  /* end index along each dimension */
      const T&                      /* value */
   );

   cmatrix<T,Syn>& subassign(
      const array<unsigned long>&,  /* start index along each dimension */
      const array<unsigned long>&,  /* step along each dimension */
      const array<unsigned long>&,  /* end index along each dimension */
      const complex<T>&             /* value */
   );

   cmatrix<T,Syn>& subassign(
      const array<unsigned long>&,  /* start index along each dimension */
      const array<unsigned long>&,  /* step along each dimension */
      const array<unsigned long>&,  /* end index along each dimension */
      const matrix<T,Syn>&          /* values */
   );
   
   cmatrix<T,Syn>& subassign(
      const array<unsigned long>&,  /* start index along each dimension */
      const array<unsigned long>&,  /* step along each dimension */
      const array<unsigned long>&,  /* end index along each dimension */
      const cmatrix<T,Syn>&         /* values */
   );
 
   cmatrix<T,Syn>& subassign(
      const array< array<unsigned long> >&,  /* indices along each dimension */
      const T&                               /* value */
   );
  
   cmatrix<T,Syn>& subassign(
      const array< array<unsigned long> >&,  /* indices along each dimension */
      const complex<T>&                      /* value */
   );

   cmatrix<T,Syn>& subassign(
      const array< array<unsigned long> >&,  /* indices along each dimension */
      const matrix<T,Syn>&                   /* values */
   );
  
   cmatrix<T,Syn>& subassign(
      const array< array<unsigned long> >&,  /* indices along each dimension */
      const cmatrix<T,Syn>&                  /* values */
   );

   /*
    * Fill with value.
    */
   cmatrix<T,Syn>& fill(const T&);
   cmatrix<T,Syn>& fill(const complex<T>&);

   /*
    * Formatted output to stream.
    */
   friend ostream& operator<< <T,Syn>(ostream&, const cmatrix<T,Syn>&);

   /*
    * Assignment operators: cmatrix-real.
    */
   cmatrix<T,Syn>& operator=(const T&);
   cmatrix<T,Syn>& operator+=(const T&);
   cmatrix<T,Syn>& operator-=(const T&);
   cmatrix<T,Syn>& operator*=(const T&);
   cmatrix<T,Syn>& operator/=(const T&);

   /*
    * Assignment operators: cmatrix-complex.
    */
   cmatrix<T,Syn>& operator=(const complex<T>&);
   cmatrix<T,Syn>& operator+=(const complex<T>&);
   cmatrix<T,Syn>& operator-=(const complex<T>&);
   cmatrix<T,Syn>& operator*=(const complex<T>&);
   cmatrix<T,Syn>& operator/=(const complex<T>&);

   /*
    * Assignment operators: cmatrix-matrix.
    */
   template <typename U, typename S> cmatrix<T,Syn>&
      operator=(const matrix<U,S>&);
   template <typename U, typename S> cmatrix<T,Syn>&
      operator+=(const matrix<U,S>&); 
   template <typename U, typename S> cmatrix<T,Syn>&
      operator-=(const matrix<U,S>&);
   template <typename U, typename S> cmatrix<T,Syn>&
      operator*=(const matrix<U,S>&);
   template <typename U, typename S> cmatrix<T,Syn>&
      operator/=(const matrix<U,S>&);

   /*
    * Assignment operators: cmatrix-cmatrix.
    */
   cmatrix<T,Syn>& operator=(const cmatrix<T,Syn>&);
   
   template <typename U, typename S> cmatrix<T,Syn>&
      operator=(const cmatrix<U,S>&);
   template <typename U, typename S> cmatrix<T,Syn>&
      operator+=(const cmatrix<U,S>&); 
   template <typename U, typename S> cmatrix<T,Syn>&
      operator-=(const cmatrix<U,S>&);
   template <typename U, typename S> cmatrix<T,Syn>&
      operator*=(const cmatrix<U,S>&);
   template <typename U, typename S> cmatrix<T,Syn>&
      operator/=(const cmatrix<U,S>&);

   /*
    * Binary cmatrix-real operators.
    */
   friend cmatrix<T,Syn> operator+ <T,Syn>(const cmatrix<T,Syn>&, const T&);
   friend cmatrix<T,Syn> operator- <T,Syn>(const cmatrix<T,Syn>&, const T&);
   friend cmatrix<T,Syn> operator* <T,Syn>(const cmatrix<T,Syn>&, const T&);
   friend cmatrix<T,Syn> operator/ <T,Syn>(const cmatrix<T,Syn>&, const T&);

   friend cmatrix<T,Syn> operator+ <T,Syn>(const T&, const cmatrix<T,Syn>&);
   friend cmatrix<T,Syn> operator- <T,Syn>(const T&, const cmatrix<T,Syn>&);
   friend cmatrix<T,Syn> operator* <T,Syn>(const T&, const cmatrix<T,Syn>&);
   friend cmatrix<T,Syn> operator/ <T,Syn>(const T&, const cmatrix<T,Syn>&);

   /*
    * Binary matrix-complex operators.
    */
   friend cmatrix<T,Syn>
      operator+ <T,Syn>(const matrix<T,Syn>&, const complex<T,Syn>&);
   friend cmatrix<T,Syn>
      operator- <T,Syn>(const matrix<T,Syn>&, const complex<T,Syn>&);
   friend cmatrix<T,Syn>
      operator* <T,Syn>(const matrix<T,Syn>&, const complex<T,Syn>&);
   friend cmatrix<T,Syn>
      operator/ <T,Syn>(const matrix<T,Syn>&, const complex<T,Syn>&);
   
   friend cmatrix<T,Syn>
      operator+ <T,Syn>(const complex<T,Syn>&, const matrix<T,Syn>&);
   friend cmatrix<T,Syn>
      operator- <T,Syn>(const complex<T,Syn>&, const matrix<T,Syn>&);
   friend cmatrix<T,Syn>
      operator* <T,Syn>(const complex<T,Syn>&, const matrix<T,Syn>&);
   friend cmatrix<T,Syn>
      operator/ <T,Syn>(const complex<T,Syn>&, const matrix<T,Syn>&);

   /*
    * Binary cmatrix-complex operators.
    */
   friend cmatrix<T,Syn>
      operator+ <T,Syn>(const cmatrix<T,Syn>&, const complex<T>&);
   friend cmatrix<T,Syn>
      operator- <T,Syn>(const cmatrix<T,Syn>&, const complex<T>&);
   friend cmatrix<T,Syn>
      operator* <T,Syn>(const cmatrix<T,Syn>&, const complex<T>&);
   friend cmatrix<T,Syn>
      operator/ <T,Syn>(const cmatrix<T,Syn>&, const complex<T>&);

   friend cmatrix<T,Syn>
      operator+ <T,Syn>(const complex<T>&, const cmatrix<T,Syn>&);
   friend cmatrix<T,Syn>
      operator- <T,Syn>(const complex<T>&, const cmatrix<T,Syn>&);
   friend cmatrix<T,Syn>
      operator* <T,Syn>(const complex<T>&, const cmatrix<T,Syn>&);
   friend cmatrix<T,Syn>
      operator/ <T,Syn>(const complex<T>&, const cmatrix<T,Syn>&);

   /*
    * Binary cmatrix-matrix operators.
    */
   friend cmatrix<T,Syn>
      operator+ <T,Syn>(const cmatrix<T,Syn>&, const matrix<T,Syn>&);
   friend cmatrix<T,Syn>
      operator- <T,Syn>(const cmatrix<T,Syn>&, const matrix<T,Syn>&);
   friend cmatrix<T,Syn>
      operator* <T,Syn>(const cmatrix<T,Syn>&, const matrix<T,Syn>&);
   friend cmatrix<T,Syn>
      operator/ <T,Syn>(const cmatrix<T,Syn>&, const matrix<T,Syn>&);

   friend cmatrix<T,Syn>
      operator+ <T,Syn>(const matrix<T,Syn>&, const cmatrix<T,Syn>&);
   friend cmatrix<T,Syn>
      operator- <T,Syn>(const matrix<T,Syn>&, const cmatrix<T,Syn>&);
   friend cmatrix<T,Syn>
      operator* <T,Syn>(const matrix<T,Syn>&, const cmatrix<T,Syn>&);
   friend cmatrix<T,Syn>
      operator/ <T,Syn>(const matrix<T,Syn>&, const cmatrix<T,Syn>&);

   /*
    * Binary cmatrix-cmatrix operators.
    */
   friend cmatrix<T,Syn>
      operator+ <T,Syn>(const cmatrix<T,Syn>&, const cmatrix<T,Syn>&);
   friend cmatrix<T,Syn>
      operator- <T,Syn>(const cmatrix<T,Syn>&, const cmatrix<T,Syn>&);
   friend cmatrix<T,Syn>
      operator* <T,Syn>(const cmatrix<T,Syn>&, const cmatrix<T,Syn>&);
   friend cmatrix<T,Syn>
      operator/ <T,Syn>(const cmatrix<T,Syn>&, const cmatrix<T,Syn>&);

   /*
    * Binary cmatrix-real comparators.
    */
   friend matrix<bool,Syn> operator== <T,Syn>(const cmatrix<T,Syn>&, const T&);
   friend matrix<bool,Syn> operator!= <T,Syn>(const cmatrix<T,Syn>&, const T&);
   friend matrix<bool,Syn> operator<  <T,Syn>(const cmatrix<T,Syn>&, const T&);
   friend matrix<bool,Syn> operator>  <T,Syn>(const cmatrix<T,Syn>&, const T&);
   friend matrix<bool,Syn> operator<= <T,Syn>(const cmatrix<T,Syn>&, const T&);
   friend matrix<bool,Syn> operator>= <T,Syn>(const cmatrix<T,Syn>&, const T&);
   
   friend matrix<bool,Syn> operator== <T,Syn>(const T&, const cmatrix<T,Syn>&);
   friend matrix<bool,Syn> operator!= <T,Syn>(const T&, const cmatrix<T,Syn>&);
   friend matrix<bool,Syn> operator<  <T,Syn>(const T&, const cmatrix<T,Syn>&);
   friend matrix<bool,Syn> operator>  <T,Syn>(const T&, const cmatrix<T,Syn>&);
   friend matrix<bool,Syn> operator<= <T,Syn>(const T&, const cmatrix<T,Syn>&);
   friend matrix<bool,Syn> operator>= <T,Syn>(const T&, const cmatrix<T,Syn>&);

   /*
    * Binary cmatrix-complex comparators.
    */
   friend matrix<bool,Syn>
      operator== <T,Syn>(const cmatrix<T,Syn>&, const complex<T>&);
   friend matrix<bool,Syn>
      operator!= <T,Syn>(const cmatrix<T,Syn>&, const complex<T>&);
   friend matrix<bool,Syn>
      operator<  <T,Syn>(const cmatrix<T,Syn>&, const complex<T>&);
   friend matrix<bool,Syn>
      operator>  <T,Syn>(const cmatrix<T,Syn>&, const complex<T>&);
   friend matrix<bool,Syn>
      operator<= <T,Syn>(const cmatrix<T,Syn>&, const complex<T>&);
   friend matrix<bool,Syn>
      operator>= <T,Syn>(const cmatrix<T,Syn>&, const complex<T>&);
   
   friend matrix<bool,Syn>
      operator== <T,Syn>(const complex<T>&, const cmatrix<T,Syn>&);
   friend matrix<bool,Syn>
      operator!= <T,Syn>(const complex<T>&, const cmatrix<T,Syn>&);
   friend matrix<bool,Syn>
      operator<  <T,Syn>(const complex<T>&, const cmatrix<T,Syn>&);
   friend matrix<bool,Syn>
      operator>  <T,Syn>(const complex<T>&, const cmatrix<T,Syn>&);
   friend matrix<bool,Syn>
      operator<= <T,Syn>(const complex<T>&, const cmatrix<T,Syn>&);
   friend matrix<bool,Syn>
      operator>= <T,Syn>(const complex<T>&, const cmatrix<T,Syn>&);

   /*
    * Binary cmatrix-matrix comparators.
    */
   friend matrix<bool,Syn>
      operator== <T,Syn>(const cmatrix<T,Syn>&, const matrix<T,Syn>&);
   friend matrix<bool,Syn>
      operator!= <T,Syn>(const cmatrix<T,Syn>&, const matrix<T,Syn>&);
   friend matrix<bool,Syn>
      operator<  <T,Syn>(const cmatrix<T,Syn>&, const matrix<T,Syn>&);
   friend matrix<bool,Syn>
      operator>  <T,Syn>(const cmatrix<T,Syn>&, const matrix<T,Syn>&);
   friend matrix<bool,Syn>
      operator<= <T,Syn>(const cmatrix<T,Syn>&, const matrix<T,Syn>&);
   friend matrix<bool,Syn>
      operator>= <T,Syn>(const cmatrix<T,Syn>&, const matrix<T,Syn>&);

   friend matrix<bool,Syn>
      operator== <T,Syn>(const matrix<T,Syn>&, const cmatrix<T,Syn>&);
   friend matrix<bool,Syn>
      operator!= <T,Syn>(const matrix<T,Syn>&, const cmatrix<T,Syn>&);
   friend matrix<bool,Syn>
      operator<  <T,Syn>(const matrix<T,Syn>&, const cmatrix<T,Syn>&);
   friend matrix<bool,Syn>
      operator>  <T,Syn>(const matrix<T,Syn>&, const cmatrix<T,Syn>&);
   friend matrix<bool,Syn>
      operator<= <T,Syn>(const matrix<T,Syn>&, const cmatrix<T,Syn>&);
   friend matrix<bool,Syn>
      operator>= <T,Syn>(const matrix<T,Syn>&, const cmatrix<T,Syn>&);

   /*
    * Binary cmatrix-cmatrix comparators.
    */
   friend matrix<bool,Syn>
      operator== <T,Syn>(const cmatrix<T,Syn>&, const cmatrix<T,Syn>&);
   friend matrix<bool,Syn>
      operator!= <T,Syn>(const cmatrix<T,Syn>&, const cmatrix<T,Syn>&);
   friend matrix<bool,Syn>
      operator<  <T,Syn>(const cmatrix<T,Syn>&, const cmatrix<T,Syn>&);
   friend matrix<bool,Syn>
      operator>  <T,Syn>(const cmatrix<T,Syn>&, const cmatrix<T,Syn>&);
   friend matrix<bool,Syn>
      operator<= <T,Syn>(const cmatrix<T,Syn>&, const cmatrix<T,Syn>&);
   friend matrix<bool,Syn>
      operator>= <T,Syn>(const cmatrix<T,Syn>&, const cmatrix<T,Syn>&);

   /*
    * Unary arithmetic operators.
    */
   cmatrix<T,Syn> operator+() const;
   cmatrix<T,Syn> operator-() const;
   cmatrix<T,Syn>& operator++();
   cmatrix<T,Syn>& operator--();
   cmatrix<T,Syn> operator++(int);
   cmatrix<T,Syn> operator--(int);
   
   /*
    * Unary element logical inversion operator.
    */
   using matrix<complex<T>,Syn>::operator!;

   /*
    * Unary element multiplicative inversion operator.
    */
   cmatrix<T,Syn> operator~() const;

   /*
    * Transcendentals on matrix elements.
    */
   friend cmatrix<T,Syn> acos<T,Syn>(const matrix<T,Syn>&);
   friend cmatrix<T,Syn> asin<T,Syn>(const matrix<T,Syn>&);
   friend cmatrix<T,Syn> acosh<T,Syn>(const matrix<T,Syn>&);
   friend cmatrix<T,Syn> atanh<T,Syn>(const matrix<T,Syn>&);

   friend cmatrix<T,Syn> sqrt<T,Syn>(const matrix<T,Syn>&);
   friend cmatrix<T,Syn> log<T,Syn>(const matrix<T,Syn>&);
   friend cmatrix<T,Syn> pow<T,Syn>(const matrix<T,Syn>&, const T&);
   friend cmatrix<T,Syn> pow<T,Syn>(const matrix<T,Syn>&, const complex<T>&);

   /*
    * Transcendentals on cmatrix elements.
    */
   friend cmatrix<T,Syn> cos<T,Syn>(const cmatrix<T,Syn>&);
   friend cmatrix<T,Syn> sin<T,Syn>(const cmatrix<T,Syn>&);
   friend cmatrix<T,Syn> tan<T,Syn>(const cmatrix<T,Syn>&);
   friend cmatrix<T,Syn> cosh<T,Syn>(const cmatrix<T,Syn>&);
   friend cmatrix<T,Syn> sinh<T,Syn>(const cmatrix<T,Syn>&);
   friend cmatrix<T,Syn> tanh<T,Syn>(const cmatrix<T,Syn>&);
   
   friend cmatrix<T,Syn> acos<T,Syn>(const cmatrix<T,Syn>&);
   friend cmatrix<T,Syn> asin<T,Syn>(const cmatrix<T,Syn>&);
   friend cmatrix<T,Syn> atan<T,Syn>(const cmatrix<T,Syn>&);
   friend cmatrix<T,Syn> acosh<T,Syn>(const cmatrix<T,Syn>&);
   friend cmatrix<T,Syn> asinh<T,Syn>(const cmatrix<T,Syn>&);
   friend cmatrix<T,Syn> atanh<T,Syn>(const cmatrix<T,Syn>&);

   friend cmatrix<T,Syn> sqrt<T,Syn>(const cmatrix<T,Syn>&);
   friend cmatrix<T,Syn> log<T,Syn>(const cmatrix<T,Syn>&);
   friend cmatrix<T,Syn> exp<T,Syn>(const cmatrix<T,Syn>&);
   friend cmatrix<T,Syn> pow<T,Syn>(const cmatrix<T,Syn>&, int);
   friend cmatrix<T,Syn> pow<T,Syn>(const cmatrix<T,Syn>&, const T&);
   friend cmatrix<T,Syn> pow<T,Syn>(const cmatrix<T,Syn>&, const complex<T>&);

   /*
    * Two-parameter arctangent of the real parts of the corresponding elements.
    */
   friend matrix<T,Syn> atan2<T,Syn>(
      const cmatrix<T,Syn>& /* y */, const cmatrix<T,Syn>& /* x */
   );

   /*
    * Magnitude of elements.
    */
   friend matrix<T,Syn> abs<T,Syn>(const cmatrix<T,Syn>&);

   /*
    * Argument of elements.
    */
   friend matrix<T,Syn> arg<T,Syn>(const cmatrix<T,Syn>&);

   /*
    * Check if matrix is real-valued.
    */
   bool is_real() const;
   
   /*
    * Real component.
    */
   matrix<T,Syn> real() const;
   
   /*
    * Imaginary component. 
    */
   matrix<T,Syn> imag() const;

   /*
    * Real and imaginary components.
    */
   void parts(matrix<T,Syn>& /* real */, matrix<T,Syn>& /* imag */) const;
   
   /*
    * Complex conjugate.
    */
   friend cmatrix<T,Syn> conj<T,Syn>(const cmatrix<T,Syn>&);
   cmatrix<T,Syn>& conj();

   /*
    * Convolution.
    */
   friend cmatrix<T,Syn>
      conv<T,Syn>(const cmatrix<T,Syn>&, const cmatrix<T,Syn>&);

   /*
    * Convolution.
    * Crop the result to be no larger than the left input matrix.
    */
   friend cmatrix<T,Syn>
      conv_crop<T,Syn>(const cmatrix<T,Syn>&, const cmatrix<T,Syn>&);

   /*
    * Convolution.
    * Crop the result to be no larger than the left input matrix and to
    * include only the central portion for which no zero padding was required.
    */
   friend cmatrix<T,Syn>
      conv_crop_strict<T,Syn>(const cmatrix<T,Syn>&, const cmatrix<T,Syn>&);

   /*
    * Dot product.
    */
   friend complex<T>
      dot<T,Syn>(const cmatrix<T,Syn>&, const cmatrix<T,Syn>&);

   /*
    * Element product.
    * Compute the product of corresponding matrix elements.
    */
   friend cmatrix<T,Syn>
      prod<T,Syn>(const cmatrix<T,Syn>&, const cmatrix<T,Syn>&);
 
   /*
    * Minimum, maximum, sum, product of all matrix elements.
    * Note that min/max operate on lexicographic ordering.
    */
   friend complex<T> min<T,Syn>(const cmatrix<T,Syn>&);
   friend complex<T> max<T,Syn>(const cmatrix<T,Syn>&);
   friend complex<T> sum<T,Syn>(const cmatrix<T,Syn>&);
   friend complex<T> prod<T,Syn>(const cmatrix<T,Syn>&);
    
   /*
    * Minimum, maximum, sum, product along specified dimension.
    * Note that min/max operate on lexicographic ordering.
    */
   friend cmatrix<T,Syn> min<T,Syn>(const cmatrix<T,Syn>&, unsigned long);
   friend cmatrix<T,Syn> max<T,Syn>(const cmatrix<T,Syn>&, unsigned long);
   friend cmatrix<T,Syn> sum<T,Syn>(const cmatrix<T,Syn>&, unsigned long);
   friend cmatrix<T,Syn> prod<T,Syn>(const cmatrix<T,Syn>&, unsigned long);
 
   /*
    * Cumulative sum, product of all matrix elements.
    */
   friend cmatrix<T,Syn> cumsum<T,Syn>(const cmatrix<T,Syn>&);
   friend cmatrix<T,Syn> cumprod<T,Syn>(const cmatrix<T,Syn>&);

   /*
    * Cumulative sum, product along specified dimension.
    */
   friend cmatrix<T,Syn> cumsum<T,Syn>(const cmatrix<T,Syn>&, unsigned long);
   friend cmatrix<T,Syn> cumprod<T,Syn>(const cmatrix<T,Syn>&, unsigned long);

   /*
    * Mean, variance of all matrix elements.
    */
   friend complex<T> mean<T,Syn>(const cmatrix<T,Syn>&);
   friend complex<T> var<T,Syn>(const cmatrix<T,Syn>&);

   /*
    * Mean, variance along specified dimension.
    */
   friend cmatrix<T,Syn> mean<T,Syn>(const cmatrix<T,Syn>&, unsigned long);
   friend cmatrix<T,Syn> var<T,Syn>(const cmatrix<T,Syn>&, unsigned long);

   /*
    * Gradient along specified dimension.
    * Optionally specify the spacing (defaults to unit spacing) of elements
    * along the dimension.
    */
   friend cmatrix<T,Syn>
      gradient<T,Syn>(const cmatrix<T,Syn>&, unsigned long);
   friend cmatrix<T,Syn>
      gradient<T,Syn>(const cmatrix<T,Syn>&, unsigned long, const complex<T>&);

   /*
    * Logical AND of elements (true iff all element(s) are nonzero).
    */
   using matrix<complex<T>,Syn>::all;
   
   /*
    * Logical OR of elements (true iff at least one element is nonzero).
    */
   using matrix<complex<T>,Syn>::some;
   
   /*
    * Logical NOR of elements (true iff all element(s) are zero).
    */
   using matrix<complex<T>,Syn>::none;
         
   /*
    * Exactly one element true (true iff exactly one element is nonzero).
    */
   using matrix<complex<T>,Syn>::one;

   /*
    * Find locations of nonzero elements.
    */
   using matrix<complex<T>,Syn>::find;
  
   /*
    * Find locations of elements with the specified value.
    */
   using matrix<complex<T>,Syn>::find_value;
   
   array<unsigned long> find_value(const T&) const;

   /*
    * Diagonal entries (returns diagonal of matrix in a vector).
    */
   friend cmatrix<T,Syn> diag<T,Syn>(const cmatrix<T,Syn>&);
   cmatrix<T,Syn> diag() const;

   /*
    * Inverse.
    * Optionally specify the numeric tolerance.
    */
   friend cmatrix<T,Syn> inv<T,Syn>(const cmatrix<T,Syn>&, const T& /* tol */);
   cmatrix<T,Syn>& inv(const T& = eps<T>() /* tol */);

   /*
    * Row echelon form.
    * Optionally specify the numeric tolerance.
    */
   friend cmatrix<T,Syn> ref<T,Syn>(const cmatrix<T,Syn>&, const T& /* tol */);
   cmatrix<T,Syn>& ref(const T& = eps<T>() /* tol */);
  
   /*
    * Reduced row echelon form.
    * Optionally specify the numeric tolerance.
    */
   friend cmatrix<T,Syn> rref<T,Syn>(const cmatrix<T,Syn>&, const T& /* tol */);
   cmatrix<T,Syn>& rref(const T& = eps<T>() /* tol */);

   /*
    * Determinant.
    * Optionally specify the numeric tolerance.
    */
   complex<T> det(const T& = eps<T>() /* tol */) const;
   
   /*
    * Rank.
    * Optionally specify the numeric tolerance.
    */
   unsigned long rank(const T& = eps<T>() /* tol */) const;
   
   /*
    * Get dimensionality of matrix.
    */
   using matrix<complex<T>,Syn>::dimensionality;

   /*
    * Get dimensions of matrix.
    */
   using matrix<complex<T>,Syn>::dimensions;

   /*
    * Transpose.
    * Reverse dimension order.
    * The input matrix is regarded as being at least two-dimensional.
    */
   friend cmatrix<T,Syn> transpose<T,Syn>(const cmatrix<T,Syn>&);
   cmatrix<T,Syn>& transpose();
 
   /*
    * Dimension shifting.
    * Shift dimensions of the matrix to the left, wrapping the leading
    * dimensions around to the right.
    */
   friend cmatrix<T,Syn> shift_dimensions<T,Syn>(
      const cmatrix<T,Syn>&,
      unsigned long                 /* number of dimensions to shift */
   );

   cmatrix<T,Syn>& shift_dimensions(
      unsigned long                 /* number of dimensions to shift */
   );
 
   /*
    * Dimension permutation.
    */
   friend cmatrix<T,Syn> permute_dimensions<T,Syn>(
      const cmatrix<T,Syn>&,
      const array<unsigned long>&   /* ordering of dimensions */
   );

   cmatrix<T,Syn>& permute_dimensions(
      const array<unsigned long>&   /* ordering of dimensions */
   );

   /*
    * Singleton dimension elimination.
    * Remove singleton dimensions, leaving the matrix elements unchanged.
    */
   friend cmatrix<T,Syn> squeeze_dimensions<T,Syn>(const cmatrix<T,Syn>&);
   cmatrix<T,Syn>& squeeze_dimensions();

   /*
    * Size.
    */
   using matrix<complex<T>,Syn>::size;
   
   /*
    * Resize (vector).
    * Return the resized vector.
    */
   friend cmatrix<T,Syn> resize<T,Syn>(
      const cmatrix<T,Syn>&,
      unsigned long                 /* length */
   );

   friend cmatrix<T,Syn> resize<T,Syn>(
      const cmatrix<T,Syn>&,
      unsigned long,                /* length */
      const T&                      /* element init value (if enlarging) */
   );
   
   friend cmatrix<T,Syn> resize<T,Syn>(
      const cmatrix<T,Syn>&,
      unsigned long,                /* length */
      const complex<T>&             /* element init value (if enlarging) */
   );
   
   /*
    * Resize (2D matrix).
    * Return the resized M x N matrix.
    */
   friend cmatrix<T,Syn> resize<T,Syn>(
      const cmatrix<T,Syn>&,
      unsigned long,                /* M */
      unsigned long                 /* N */
   );
   
   friend cmatrix<T,Syn> resize<T,Syn>(
      const cmatrix<T,Syn>&,
      unsigned long,                /* M */
      unsigned long,                /* N */
      const T&                      /* element init value (if enlarging) */
   );
   
   friend cmatrix<T,Syn> resize<T,Syn>(
      const cmatrix<T,Syn>&,
      unsigned long,                /* M */
      unsigned long,                /* N */
      const complex<T>&             /* element init value (if enlarging) */
   );
   
   /*
    * Resize (multi-dimensional matrix).
    * Return the resized matrix.
    * The dimensionality of the matrix must not change.
    */
   friend cmatrix<T,Syn> resize<T,Syn>(
      const cmatrix<T,Syn>&,
      const array<unsigned long>&   /* dimensions */
   );
   
   friend cmatrix<T,Syn> resize<T,Syn>(
      const cmatrix<T,Syn>&,
      const array<unsigned long>&,  /* dimensions */
      const T&                      /* element init value (if enlarging) */
   );
   
   friend cmatrix<T,Syn> resize<T,Syn>(
      const cmatrix<T,Syn>&,
      const array<unsigned long>&,  /* dimensions */
      const complex<T>&             /* element init value (if enlarging) */
   );

   /*
    * Resize (vector).
    * The original vector is resized to the specified length.
    */
   cmatrix<T,Syn>& resize(
      unsigned long                 /* length */
   );

   cmatrix<T,Syn>& resize(
      unsigned long,                /* length */
      const T&                      /* element init value (if enlarging) */
   );
   
   cmatrix<T,Syn>& resize(
      unsigned long,                /* length */
      const complex<T>&             /* element init value (if enlarging) */
   );

   /*
    * Resize (2D matrix).
    * The original 2D matrix is resized to be M x N.
    */
   cmatrix<T,Syn>& resize(
      unsigned long,                /* M */
      unsigned long                 /* N */
   );
   
   cmatrix<T,Syn>& resize(
      unsigned long,                /* M */
      unsigned long,                /* N */
      const T&                      /* element init value (if enlarging) */
   );
   
   cmatrix<T,Syn>& resize(
      unsigned long,                /* M */
      unsigned long,                /* N */
      const complex<T>&             /* element init value (if enlarging) */
   );
   
   /*
    * Resize (multi-dimensional matrix).
    * The original matrix is resized to the given dimensions.
    * The dimensionality of the matrix must not change.
    */
   cmatrix<T,Syn>& resize(
      const array<unsigned long>&   /* dimensions */
   );
   
   cmatrix<T,Syn>& resize(
      const array<unsigned long>&,  /* dimensions */
      const T&                      /* element init value (if enlarging) */
   );
   
   cmatrix<T,Syn>& resize(
      const array<unsigned long>&,  /* dimensions */
      const complex<T>&             /* element init value (if enlarging) */
   );
   
   /*
    * Reshape into a vector.
    */
   friend cmatrix<T,Syn> vector<T,Syn>(const cmatrix<T,Syn>&);
   
   /*
    * Reshape into 2D M x N matrix.
    */
   friend cmatrix<T,Syn> reshape<T,Syn>(
      const cmatrix<T,Syn>&,
      unsigned long,                /* M */
      unsigned long                 /* N */
   );
   
   /*
    * Reshape into multi-dimensional matrix.
    */
   friend cmatrix<T,Syn> reshape<T,Syn>(
      const cmatrix<T,Syn>&,
      const array<unsigned long>&   /* dimensions */
   );
   
   /*
    * Reshape (in place) into a vector.
    */
   cmatrix<T,Syn>& vector();
  
   /*
    * Reshape (in place) into 2D M x N matrix.
    */
   cmatrix<T,Syn>& reshape(
      unsigned long,                /* M */
      unsigned long                 /* N */
   );
   
   /*
    * Reshape (in place) into multi-dimensional matrix.
    */
   cmatrix<T,Syn>& reshape(
      const array<unsigned long>&   /* dimensions */
   );
   
   /*
    * Replicate matrix.
    */
   friend cmatrix<T,Syn> repmat<T,Syn>(
      const cmatrix<T,Syn>&,
      unsigned long,                /* M */
      unsigned long                 /* N */
   );
      
   friend cmatrix<T,Syn> repmat<T,Syn>(
      const cmatrix<T,Syn>&, 
      const array<unsigned long>&   /* dimensions */
   );
   
   /*
    * Concatenate matrices.
    * Throw an exception (ex_matrix_dimension_mismatch) if the matrices
    * cannot be concatenated along the specified dimension.
    */
   friend cmatrix<T,Syn> vertcat<T,Syn>(
      const cmatrix<T,Syn>&, const cmatrix<T,Syn>&
   );
   friend cmatrix<T,Syn> horzcat<T,Syn>(
      const cmatrix<T,Syn>&, const cmatrix<T,Syn>&
   );
   friend cmatrix<T,Syn> concat<T,Syn>(
      const cmatrix<T,Syn>&, 
      const cmatrix<T,Syn>&, 
      unsigned long           /* dimension along which to concatenate */
   );
   
   /*
    * Reverse element order (treating the matrix as a linear array).
    */
   cmatrix<T,Syn>& reverse();

   /*
    * Reverse element order along the specified dimension.
    */
   cmatrix<T,Syn>& reverse(unsigned long /* dimension */);

   /*
    * Random permutation.
    */
   using matrix<complex<T>,Syn>::randperm;

   /*
    * Sorting.
    */
   using matrix<complex<T>,Syn>::sort;
   using matrix<complex<T>,Syn>::sort_idx;

   /*
    * Uniqueness.
    */
   using matrix<complex<T>,Syn>::unique;
   using matrix<complex<T>,Syn>::unique_idx;

protected:
   /************************************************************************
    * Matrix helper functions/methods.
    * Note: These functions do not lock their matrix arguments.
    ************************************************************************/

   /*
    * Convert a matrix<complex<T>,Syn> object into a cmatrix<T,Syn> object.
    * The original matrix<complex<T>,Syn> object is destroyed in the process.
    */
   static cmatrix<T,Syn> to_cmatrix(matrix<complex<T>,Syn>&);

   /*
    * Output classname to stream.
    */
   virtual ostream& output_classname(ostream&) const;

   /*
    * Output matrix element at given index to stream.
    */
   virtual ostream& output_element(ostream&, unsigned long) const;
};

/***************************************************************************
 * Constructors and destructor.
 ***************************************************************************/

/*
 * Constructor.
 * Return the empty (0 x 0) matrix.
 */
template <typename T, typename Syn>
cmatrix<T,Syn>::cmatrix()
 : matrix<complex<T>,Syn>()
{ }

/*
 * Constructors for single element (1 x 1) matrix.
 */
template <typename T, typename Syn>
cmatrix<T,Syn>::cmatrix(const T& t)
 : matrix<complex<T>,Syn>(complex<T>(t))
{ }

template <typename T, typename Syn>
cmatrix<T,Syn>::cmatrix(const complex<T>& z)
 : matrix<complex<T>,Syn>(z)
{ }

/*
 * Constructors for two-dimensional M x N matrices.
 */
template <typename T, typename Syn>
cmatrix<T,Syn>::cmatrix(unsigned long M, unsigned long N)
 : matrix<complex<T>,Syn>(M, N)
{ }

template <typename T, typename Syn>
cmatrix<T,Syn>::cmatrix(unsigned long M, unsigned long N, const T& t)
 : matrix<complex<T>,Syn>(M, N, complex<T>(t))
{ }

template <typename T, typename Syn>
cmatrix<T,Syn>::cmatrix(unsigned long M, unsigned long N, const complex<T>& z)
 : matrix<complex<T>,Syn>(M, N, z)
{ }

/*
 * Constructors for multi-dimensional matrices.
 */
template <typename T, typename Syn>
cmatrix<T,Syn>::cmatrix(const array<unsigned long>& dims)
 : matrix<complex<T>,Syn>(dims)
{ }

template <typename T, typename Syn>
cmatrix<T,Syn>::cmatrix(const array<unsigned long>& dims, const T& t)
 : matrix<complex<T>,Syn>(dims, complex<T>(t))
{ }

template <typename T, typename Syn>
cmatrix<T,Syn>::cmatrix(const array<unsigned long>& dims, const complex<T>& z)
 : matrix<complex<T>,Syn>(dims, z)
{ }

/*
 * Constructor from real-valued matrix.
 */
template <typename T, typename Syn>
cmatrix<T,Syn>::cmatrix(const matrix<T,Syn>& m)
 : matrix<complex<T>,Syn>(m)
{ }

/*
 * Constructor from real and imaginary components.
 */
template <typename T, typename Syn>
cmatrix<T,Syn>::cmatrix(
   const matrix<T,Syn>& m_real, const matrix<T,Syn>& m_imag)
 : matrix<complex<T>,Syn>()
{
   auto_read_read_lock<const Syn> rrlock(m_real, m_imag);
   matrix<T,Syn>::assert_dims_equal(m_real._dims, m_imag._dims);
   cmatrix<T,Syn> m(m_real._dims);
   for (unsigned long n = 0; n < this->_size; n++)
      m._data[n] = complex<T>(m_real._data[n], m_imag._data[n]);
   matrix<complex<T>,Syn>::swap(*this, m);
}

/*
 * Copy constructors.
 */
template <typename T, typename Syn>
cmatrix<T,Syn>::cmatrix(const cmatrix<T,Syn>& m)
 : matrix<complex<T>,Syn>(static_cast<const matrix<complex<T>,Syn>&>(m))
{ }

template <typename T, typename Syn>
template <typename U, typename S>
cmatrix<T,Syn>::cmatrix(const cmatrix<U,S>& m)
 : matrix<complex<T>,Syn>(static_cast<const matrix<complex<U>,S>&>(m))
{ }

/*
 * Destructor.
 */
template <typename T, typename Syn>
cmatrix<T,Syn>::~cmatrix() {
   /* do nothing */
}

/***************************************************************************
 * Named constructors.
 ***************************************************************************/

/*
 * Empty matrix (0 x 0).
 */
template <typename T, typename Syn>
cmatrix<T,Syn> cmatrix<T,Syn>::empty() {
   return cmatrix<T,Syn>();
}

/*
 * Identity matrix.
 * Create and return an N x N identity matrix.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> cmatrix<T,Syn>::eye(unsigned long N) {
   cmatrix<T,Syn> m(N, N);
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
cmatrix<T,Syn> cmatrix<T,Syn>::zeros(unsigned long M, unsigned long N) {
   return cmatrix<T,Syn>(M, N);
}

/*
 * Zeros matrix.
 * Create and return a matrix of zeros with the given dimensions.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> cmatrix<T,Syn>::zeros(const array<unsigned long>& dims) {
   return cmatrix<T,Syn>(dims);
}

/*
 * Ones matrix.
 * Create and return an M x N matrix of ones.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> cmatrix<T,Syn>::ones(unsigned long M, unsigned long N) {
   return cmatrix<T,Syn>(M, N, T(1));
}

/*
 * Ones matrix.
 * Create and return a matrix of ones with the given dimensions.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> cmatrix<T,Syn>::ones(const array<unsigned long>& dims) {
   return cmatrix<T,Syn>(dims, T(1));
}

/*
 * Random matrix.
 * Create and return an M x N matrix of random numbers drawn from the
 * given random number generator.  If no generator is specified, then
 * the matrix elements are uniformly distributed random numbers in the 
 * interval [0, 1] (inclusive).
 */ 
template <typename T, typename Syn>
cmatrix<T,Syn> cmatrix<T,Syn>::random(
   unsigned long M, unsigned long N)
{
   rand_gen_uniform<T> r;
   return cmatrix<T,Syn>::random(M, N, r);
}

template <typename T, typename Syn>
cmatrix<T,Syn> cmatrix<T,Syn>::random(
   unsigned long M, unsigned long N, rand_gen<T>& r)
{
   cmatrix<T,Syn> m(M, N);
   for (unsigned long n = 0; n < m._size; n++)
      m._data[n] = r.generate();
   return m;
}

template <typename T, typename Syn>
cmatrix<T,Syn> cmatrix<T,Syn>::random(
   unsigned long M, unsigned long N, rand_gen< complex<T> >& r)
{
   cmatrix<T,Syn> m(M, N);
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
cmatrix<T,Syn> cmatrix<T,Syn>::random(
   const array<unsigned long>& dims)
{
   rand_gen_uniform<T> r;
   return cmatrix<T,Syn>::random(dims, r);
}

template <typename T, typename Syn>
cmatrix<T,Syn> cmatrix<T,Syn>::random(
   const array<unsigned long>& dims, rand_gen<T>& r)
{
   cmatrix<T,Syn> m(dims);
   for (unsigned long n = 0; n < m._size; n++)
      m._data[n] = r.generate();
   return m;
}

template <typename T, typename Syn>
cmatrix<T,Syn> cmatrix<T,Syn>::random(
   const array<unsigned long>& dims, rand_gen< complex<T> >& r)
{
   cmatrix<T,Syn> m(dims);
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
cmatrix<T,Syn> cmatrix<T,Syn>::ramp(
   const T& start, const T& end)
{
   const T step(1);
   return cmatrix<T,Syn>::ramp(start, step, end);
}

template <typename T, typename Syn>
cmatrix<T,Syn> cmatrix<T,Syn>::ramp(
   const T& start, const T& step, const T& end)
{
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
   cmatrix<T,Syn> m(n_els, 1);
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
cmatrix<T,Syn> cmatrix<T,Syn>::diagonal(
   const matrix<T,Syn>& dgnl, unsigned long d)
{
   matrix<complex<T>,Syn> m = matrix<complex<T>,Syn>::diagonal(
      matrix<complex<T>,Syn>(dgnl), d
   );
   return cmatrix<T,Syn>::to_cmatrix(m);
}

template <typename T, typename Syn>
cmatrix<T,Syn> cmatrix<T,Syn>::diagonal(
   const cmatrix<T,Syn>& dgnl, unsigned long d)
{
   matrix<complex<T>,Syn> m = matrix<complex<T>,Syn>::diagonal(
      dgnl, d
   );
   return cmatrix<T,Syn>::to_cmatrix(m);
}

/***************************************************************************
 * Conversion to cmatrix.
 ***************************************************************************/

/*
 * Convert a matrix<complex<T>,Syn> object into a cmatrix<T,Syn> object.
 * The original matrix<complex<T>,Syn> object is destroyed in the process.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> cmatrix<T,Syn>::to_cmatrix(matrix<complex<T>,Syn>& m) {
   cmatrix<T,Syn> cm;
   matrix<complex<T>,Syn>::swap(cm, m);
   return cm;
}

/***************************************************************************
 * Serialization.
 ***************************************************************************/

/*
 * Serialize.
 */
template <typename T, typename Syn>
void cmatrix<T,Syn>::serialize(
   serial_output_stream& s, const serializer<T>& slzr) const
{
   const container_serializer<complex<T>,T> cmplx_slzr(slzr);
   this->matrix<complex<T>,Syn>::serialize(s, cmplx_slzr);
}

/*
 * Deserialize.
 */
template <typename T, typename Syn>
auto_ptr< cmatrix<T,Syn> > cmatrix<T,Syn>::deserialize(
   serial_input_stream& s, const serializer<T>& slzr)
{
   const container_serializer<complex<T>,T> cmplx_slzr(slzr);
   auto_ptr< cmatrix<T,Syn> > cm(new cmatrix<T,Syn>());
   auto_ptr< matrix<complex<T>,Syn> > m =
      matrix<complex<T>,Syn>::deserialize(s, cmplx_slzr);
   matrix<complex<T>,Syn>::swap(*cm, *m);
   return cm;
}

/***************************************************************************
 * Subarray.
 ***************************************************************************/

/*
 * Subarray.
 * Return vector of elements in the specified linear index range.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> cmatrix<T,Syn>::subarray(
   unsigned long start,
   unsigned long end) const
{
   matrix<complex<T>,Syn> m = this->matrix<complex<T>,Syn>::subarray(
      start, end
   );
   return cmatrix<T,Syn>::to_cmatrix(m);
}

template <typename T, typename Syn>
cmatrix<T,Syn> cmatrix<T,Syn>::subarray(
   unsigned long start, 
   unsigned long step,
   unsigned long end) const
{
   matrix<complex<T>,Syn> m = this->matrix<complex<T>,Syn>::subarray(
      start, step, end
   );
   return cmatrix<T,Syn>::to_cmatrix(m);
}

template <typename T, typename Syn>
cmatrix<T,Syn> cmatrix<T,Syn>::subarray(
   const array<unsigned long>& inds) const
{
   matrix<complex<T>,Syn> m = this->matrix<complex<T>,Syn>::subarray(inds);
   return cmatrix<T,Syn>::to_cmatrix(m);
}

/***************************************************************************
 * Submatrix.
 ***************************************************************************/
 
/*
 * Submatrix.
 * Return matrix of elements in the specified range.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> cmatrix<T,Syn>::submatrix(
   const array<unsigned long>& start,
   const array<unsigned long>& end) const
{
   matrix<complex<T>,Syn> m = this->matrix<complex<T>,Syn>::submatrix(
      start, end
   );
   return cmatrix<T,Syn>::to_cmatrix(m);
}

template <typename T, typename Syn>
cmatrix<T,Syn> cmatrix<T,Syn>::submatrix(
   const array<unsigned long>& start,
   const array<unsigned long>& step,
   const array<unsigned long>& end) const
{
   matrix<complex<T>,Syn> m = this->matrix<complex<T>,Syn>::submatrix(
      start, step, end
   );
   return cmatrix<T,Syn>::to_cmatrix(m);
}

template <typename T, typename Syn>
cmatrix<T,Syn> cmatrix<T,Syn>::submatrix(
   const array< array<unsigned long> >& range) const
{
   matrix<complex<T>,Syn> m = this->matrix<complex<T>,Syn>::submatrix(range);
   return cmatrix<T,Syn>::to_cmatrix(m);
}

/***************************************************************************
 * Subassign.
 ***************************************************************************/

/*
 * Subarray assignment.
 * Assign value(s) to elements in the specified linear index range.
 */
template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::subassign(
   unsigned long start,
   unsigned long end,
   const T& t)
{
   return this->subassign(start, end, complex<T>(t));
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::subassign(
   unsigned long start,
   unsigned long end,
   const complex<T>& z)
{
   this->matrix<complex<T>,Syn>::subassign(start, end, z);
   return *this;
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::subassign(
   unsigned long start,
   unsigned long end,
   const matrix<T,Syn>& values)
{
   return this->subassign(start, end, cmatrix<T,Syn>(values));
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::subassign(
   unsigned long start,
   unsigned long end,
   const cmatrix<T,Syn>& values)
{
   this->matrix<complex<T>,Syn>::subassign(
      start, end, static_cast<const matrix<complex<T>,Syn>&>(values)
   );
   return *this;
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::subassign(
   unsigned long start,
   unsigned long step,
   unsigned long end,
   const T& t)
{
   return this->subassign(start, step, end, complex<T>(t));
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::subassign(
   unsigned long start,
   unsigned long step,
   unsigned long end,
   const complex<T>& z)
{
   this->matrix<complex<T>,Syn>::subassign(start, step, end, z);
   return *this;
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::subassign(
   unsigned long start,
   unsigned long step,
   unsigned long end,
   const matrix<T,Syn>& values)
{
   return this->subassign(start, step, end, cmatrix<T,Syn>(values));
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::subassign(
   unsigned long start,
   unsigned long step,
   unsigned long end,
   const cmatrix<T,Syn>& values)
{
   this->matrix<complex<T>,Syn>::subassign(
      start, step, end, static_cast<const matrix<complex<T>,Syn>&>(values)
   );
   return *this;
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::subassign(
   const array<unsigned long>& inds,
   const T& t)
{
   this->subassign(inds, complex<T>(t));
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::subassign(
   const array<unsigned long>& inds,
   const complex<T>& z)
{
   this->matrix<complex<T>,Syn>::subassign(inds, z);
   return *this;
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::subassign(
   const array<unsigned long>& inds,
   const matrix<T,Syn>& values)
{
   return this->subassign(inds, cmatrix<T,Syn>(values));
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::subassign(
   const array<unsigned long>& inds,
   const cmatrix<T,Syn>& values)
{
   this->matrix<complex<T>,Syn>::subassign(
      inds, static_cast<const matrix<complex<T>,Syn>&>(values)
   );
   return *this;
}

/*
 * Submatrix assignment.
 * Assign value(s) to elements in the specified range.
 */
template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::subassign(
   const array<unsigned long>& start,
   const array<unsigned long>& end,
   const T& t)
{
   return this->subassign(start, end, complex<T>(t));
}
   
template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::subassign(
   const array<unsigned long>& start,
   const array<unsigned long>& end,
   const complex<T>& z)
{
   this->matrix<complex<T>,Syn>::subassign(start, end, z);
   return *this;
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::subassign(
   const array<unsigned long>& start,
   const array<unsigned long>& end,
   const matrix<T,Syn>& values)
{
   return this->subassign(start, end, cmatrix<T,Syn>(values));
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::subassign(
   const array<unsigned long>& start,
   const array<unsigned long>& end,
   const cmatrix<T,Syn>& values)
{
   this->matrix<complex<T>,Syn>::subassign(
      start, end, static_cast<const matrix<complex<T>,Syn>&>(values)
   );
   return *this;
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::subassign(
   const array<unsigned long>& start,
   const array<unsigned long>& step,
   const array<unsigned long>& end,
   const T& t)
{
   return this->subassign(start, step, end, complex<T>(t));
}
   
template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::subassign(
   const array<unsigned long>& start,
   const array<unsigned long>& step,
   const array<unsigned long>& end,
   const complex<T>& z)
{
   this->matrix<complex<T>,Syn>::subassign(start, step, end, z);
   return *this;
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::subassign(
   const array<unsigned long>& start,
   const array<unsigned long>& step,
   const array<unsigned long>& end,
   const matrix<T,Syn>& values)
{
   return this->subassign(start, step, end, cmatrix<T,Syn>(values));
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::subassign(
   const array<unsigned long>& start,
   const array<unsigned long>& step,
   const array<unsigned long>& end,
   const cmatrix<T,Syn>& values)
{
   this->matrix<complex<T>,Syn>::subassign(
      start, step, end, static_cast<const matrix<complex<T>,Syn>&>(values)
   );
   return *this;
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::subassign(
   const array< array<unsigned long> >& range,
   const T& t)
{
   return this->subassign(range, complex<T>(t));
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::subassign(
   const array< array<unsigned long> >& range,
   const complex<T>& z)
{
   this->matrix<complex<T>,Syn>::subassign(range, z);
   return *this;
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::subassign(
   const array< array<unsigned long> >& range,
   const matrix<T,Syn>& values)
{
   return this->subassign(range, cmatrix<T,Syn>(values));
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::subassign(
   const array< array<unsigned long> >& range,
   const cmatrix<T,Syn>& values)
{
   this->matrix<complex<T>,Syn>::subassign(
      range, static_cast<const matrix<complex<T>,Syn>&>(values)
   );
   return *this;
}

/***************************************************************************
 * Fill.
 ***************************************************************************/

/*
 * Fill with value.
 */
template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::fill(const T& t) {
   return this->fill(complex<T>(t));
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::fill(const complex<T>& z) {
   this->matrix<complex<T>,Syn>::fill(z);
   return *this;
}

/***************************************************************************
 * I/O.
 ***************************************************************************/

/*
 * Formatted output to stream.
 */
template <typename T, typename Syn>
ostream& operator<<(ostream& os, const cmatrix<T,Syn>& m) {
   return (os << static_cast<const matrix<complex<T>,Syn>&>(m));
}

/*
 * Output classname to stream.
 */
template <typename T, typename Syn>
ostream& cmatrix<T,Syn>::output_classname(ostream& os) const {
   os << "cmatrix";
   return os;
}

/*
 * Output matrix element at given index to stream.
 */
template <typename T, typename Syn>
ostream& cmatrix<T,Syn>::output_element(ostream& os, unsigned long n) const {
   /* get real and imaginary components */
   T t_real = this->_data[n].real();
   T t_imag = this->_data[n].imag();
   /* output real part to string */
   ostringstream s0;
   s0 << io::streams::iomanip::setw(10)
      << io::streams::iomanip::setprecision(4)
      << io::streams::iomanip::setiosflags(ios::right)
      << t_real;
   /* output imaginary part to string */
   if (t_imag != T()) {
      s0 << ((t_imag < T()) ? " - " : " + ")
         << io::streams::iomanip::setprecision(4)
         << io::streams::iomanip::resetiosflags(ios::right)
         << io::streams::iomanip::setiosflags(ios::left)
         << abs(t_imag) << 'i';
   }
   /* create string containing element in fixed field size */
   ostringstream s1;
   s1 << io::streams::iomanip::setw(24)
      << io::streams::iomanip::setiosflags(ios::left)
      << s0.str();
   /* output string */
   os << s1.str();
   return os;
}

/***************************************************************************
 * Assignment operators: cmatrix-real.
 ***************************************************************************/

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::operator=(const T& t) {
   return (*this = complex<T>(t));
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::operator+=(const T& t) {
   return (*this += complex<T>(t));
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::operator-=(const T& t) {
   return (*this -= complex<T>(t));
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::operator*=(const T& t) {
   return (*this *= complex<T>(t));
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::operator/=(const T& t) {
   return (*this /= complex<T>(t));
}

/***************************************************************************
 * Assignment operators: cmatrix-complex.
 ***************************************************************************/

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::operator=(const complex<T>& z) {
   this->matrix<complex<T>,Syn>::operator=(z);
   return *this;
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::operator+=(const complex<T>& z) {
   this->matrix<complex<T>,Syn>::operator+=(z);
   return *this;
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::operator-=(const complex<T>& z) {
   this->matrix<complex<T>,Syn>::operator-=(z);
   return *this;
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::operator*=(const complex<T>& z) {
   this->matrix<complex<T>,Syn>::operator*=(z);
   return *this;
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::operator/=(const complex<T>& z) {
   this->matrix<complex<T>,Syn>::operator/=(z);
   return *this;
}

/***************************************************************************
 * Assignment operators: cmatrix-matrix.
 ***************************************************************************/

template <typename T, typename Syn>
template <typename U, typename S>
cmatrix<T,Syn>& cmatrix<T,Syn>::operator=(const matrix<U,S>& m) {
   this->matrix<complex<T>,Syn>::operator=(m);
   return *this;
}

template <typename T, typename Syn>
template <typename U, typename S>
cmatrix<T,Syn>& cmatrix<T,Syn>::operator+=(const matrix<U,S>& m) {
   this->matrix<complex<T>,Syn>::operator+=(m);
   return *this;
}

template <typename T, typename Syn>
template <typename U, typename S>
cmatrix<T,Syn>& cmatrix<T,Syn>::operator-=(const matrix<U,S>& m) {
   this->matrix<complex<T>,Syn>::operator-=(m);
   return *this;
}

template <typename T, typename Syn>
template <typename U, typename S>
cmatrix<T,Syn>& cmatrix<T,Syn>::operator*=(const matrix<U,S>& m) {
   this->matrix<complex<T>,Syn>::operator*=(m);
   return *this;
}

template <typename T, typename Syn>
template <typename U, typename S>
cmatrix<T,Syn>& cmatrix<T,Syn>::operator/=(const matrix<U,S>& m) {
   this->matrix<complex<T>,Syn>::operator/=(m);
   return *this;
}

/***************************************************************************
 * Assignment operators: cmatrix-cmatrix.
 ***************************************************************************/

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::operator=(const cmatrix<T,Syn>& m) {
   this->matrix<complex<T>,Syn>::operator=(
      static_cast<const matrix<complex<T>,Syn>&>(m)
   );
   return *this;
}

template <typename T, typename Syn>
template <typename U, typename S>
cmatrix<T,Syn>& cmatrix<T,Syn>::operator=(const cmatrix<U,S>& m) {
   this->matrix<complex<T>,Syn>::operator=(
      static_cast<const matrix<complex<U>,S>&>(m)
   );
   return *this;
}

template <typename T, typename Syn>
template <typename U, typename S>
cmatrix<T,Syn>& cmatrix<T,Syn>::operator+=(const cmatrix<U,S>& m) {
   this->matrix<complex<T>,Syn>::operator+=(
      static_cast<const matrix<complex<U>,S>&>(m)
   );
   return *this;
}

template <typename T, typename Syn>
template <typename U, typename S>
cmatrix<T,Syn>& cmatrix<T,Syn>::operator-=(const cmatrix<U,S>& m) {
   this->matrix<complex<T>,Syn>::operator-=(
      static_cast<const matrix<complex<U>,S>&>(m)
   );
   return *this;
}

template <typename T, typename Syn>
template <typename U, typename S>
cmatrix<T,Syn>& cmatrix<T,Syn>::operator*=(const cmatrix<U,S>& m) {
   this->matrix<complex<T>,Syn>::operator*=(
      static_cast<const matrix<complex<U>,S>&>(m)
   );
   return *this;
}

template <typename T, typename Syn>
template <typename U, typename S>
cmatrix<T,Syn>& cmatrix<T,Syn>::operator/=(const cmatrix<U,S>& m) {
   this->matrix<complex<T>,Syn>::operator/=(
      static_cast<const matrix<complex<U>,S>&>(m)
   );
   return *this;
}

/***************************************************************************
 * Binary cmatrix-real operators.
 ***************************************************************************/
 
template <typename T, typename Syn>
cmatrix<T,Syn> operator+(const cmatrix<T,Syn>& m, const T& t) {
   return (cmatrix<T,Syn>(m) += t);
}

template <typename T, typename Syn>
cmatrix<T,Syn> operator-(const cmatrix<T,Syn>& m, const T& t) {
   return (cmatrix<T,Syn>(m) -= t);
}

template <typename T, typename Syn>
cmatrix<T,Syn> operator*(const cmatrix<T,Syn>& m, const T& t) {
   return (cmatrix<T,Syn>(m) *= t);
}

template <typename T, typename Syn>
cmatrix<T,Syn> operator/(const cmatrix<T,Syn>& m, const T& t) {
   return (cmatrix<T,Syn>(m) /= t);
}

template <typename T, typename Syn>
cmatrix<T,Syn> operator+(const T& t, const cmatrix<T,Syn>& m) {
   return (cmatrix<T,Syn>(m) += t);
}

template <typename T, typename Syn>
cmatrix<T,Syn> operator-(const T& t, const cmatrix<T,Syn>& m) {
   cmatrix<T,Syn> m_result = (-m);
   m_result += t;
   return m_result;
}

template <typename T, typename Syn>
cmatrix<T,Syn> operator*(const T& t, const cmatrix<T,Syn>& m) {
   return (cmatrix<T,Syn>(m) *= t);
}

template <typename T, typename Syn>
cmatrix<T,Syn> operator/(const T& t, const cmatrix<T,Syn>& m) {
   cmatrix<T,Syn> m_result = inv(m);
   m_result *= t;
   return m_result;
}

/***************************************************************************
 * Binary matrix-complex operators.
 ***************************************************************************/

template <typename T, typename Syn>
cmatrix<T,Syn> operator+(const matrix<T,Syn>& m, const complex<T,Syn>& z) {
   return (cmatrix<T,Syn>(m) += z);
}

template <typename T, typename Syn>
cmatrix<T,Syn> operator-(const matrix<T,Syn>& m, const complex<T,Syn>& z) {
   return (cmatrix<T,Syn>(m) -= z);
}

template <typename T, typename Syn>
cmatrix<T,Syn> operator*(const matrix<T,Syn>& m, const complex<T,Syn>& z) {
   return (cmatrix<T,Syn>(m) *= z);
}

template <typename T, typename Syn>
cmatrix<T,Syn> operator/(const matrix<T,Syn>& m, const complex<T,Syn>& z) {
   return (cmatrix<T,Syn>(m) /= z);
}

template <typename T, typename Syn>
cmatrix<T,Syn> operator+(const complex<T,Syn>& z, const matrix<T,Syn>& m) {
   return (cmatrix<T,Syn>(m) += z);
}

template <typename T, typename Syn>
cmatrix<T,Syn> operator-(const complex<T,Syn>& z, const matrix<T,Syn>& m) {
   cmatrix<T,Syn> m_result(-m);
   m_result += z;
   return m_result;
}

template <typename T, typename Syn>
cmatrix<T,Syn> operator*(const complex<T,Syn>& z, const matrix<T,Syn>& m) {
   return (cmatrix<T,Syn>(m) *= z);
}

template <typename T, typename Syn>
cmatrix<T,Syn> operator/(const complex<T,Syn>& z, const matrix<T,Syn>& m) {
   cmatrix<T,Syn> m_result(inv(m));
   m_result *= z;
   return m_result;
}

/***************************************************************************
 * Binary cmatrix-complex operators.
 ***************************************************************************/

template <typename T, typename Syn>
cmatrix<T,Syn> operator+(const cmatrix<T,Syn>& m, const complex<T>& z) {
   return (cmatrix<T,Syn>(m) += z);
}

template <typename T, typename Syn>
cmatrix<T,Syn> operator-(const cmatrix<T,Syn>& m, const complex<T>& z) {
   return (cmatrix<T,Syn>(m) -= z);
}

template <typename T, typename Syn>
cmatrix<T,Syn> operator*(const cmatrix<T,Syn>& m, const complex<T>& z) {
   return (cmatrix<T,Syn>(m) *= z);
}

template <typename T, typename Syn>
cmatrix<T,Syn> operator/(const cmatrix<T,Syn>& m, const complex<T>& z) {
   return (cmatrix<T,Syn>(m) /= z);
}

template <typename T, typename Syn>
cmatrix<T,Syn> operator+(const complex<T>& z, const cmatrix<T,Syn>& m) {
   return (cmatrix<T,Syn>(m) += z);
}

template <typename T, typename Syn>
cmatrix<T,Syn> operator-(const complex<T>& z, const cmatrix<T,Syn>& m) {
   cmatrix<T,Syn> m_result = (-m);
   m_result += z;
   return m_result;
}

template <typename T, typename Syn>
cmatrix<T,Syn> operator*(const complex<T>& z, const cmatrix<T,Syn>& m) {
   return (cmatrix<T,Syn>(m) *= z);
}

template <typename T, typename Syn>
cmatrix<T,Syn> operator/(const complex<T>& z, const cmatrix<T,Syn>& m) {
   cmatrix<T,Syn> m_result = inv(m);
   m_result *= z;
   return m_result;
}

/***************************************************************************
 * Binary cmatrix-matrix operators.
 ***************************************************************************/

template <typename T, typename Syn>
cmatrix<T,Syn> operator+(const cmatrix<T,Syn>& m0, const matrix<T,Syn>& m1) {
   return (cmatrix<T,Syn>(m0) += m1);
}

template <typename T, typename Syn>
cmatrix<T,Syn> operator-(const cmatrix<T,Syn>& m0, const matrix<T,Syn>& m1) {
   return (cmatrix<T,Syn>(m0) -= m1);
}

template <typename T, typename Syn>
cmatrix<T,Syn> operator*(const cmatrix<T,Syn>& m0, const matrix<T,Syn>& m1) {
   return (cmatrix<T,Syn>(m0) *= m1);
}

template <typename T, typename Syn>
cmatrix<T,Syn> operator/(const cmatrix<T,Syn>& m0, const matrix<T,Syn>& m1) {
   return (cmatrix<T,Syn>(m0) /= m1);
}

template <typename T, typename Syn>
cmatrix<T,Syn> operator+(const matrix<T,Syn>& m0, const cmatrix<T,Syn>& m1) {
   return (cmatrix<T,Syn>(m0) += m1);
}

template <typename T, typename Syn>
cmatrix<T,Syn> operator-(const matrix<T,Syn>& m0, const cmatrix<T,Syn>& m1) {
   return (cmatrix<T,Syn>(m0) -= m1);
}

template <typename T, typename Syn>
cmatrix<T,Syn> operator*(const matrix<T,Syn>& m0, const cmatrix<T,Syn>& m1) {
   return (cmatrix<T,Syn>(m0) *= m1);
}

template <typename T, typename Syn>
cmatrix<T,Syn> operator/(const matrix<T,Syn>& m0, const cmatrix<T,Syn>& m1) {
   return (cmatrix<T,Syn>(m0) /= m1);
}

/***************************************************************************
 * Binary cmatrix-cmatrix operators.
 ***************************************************************************/

template <typename T, typename Syn>
cmatrix<T,Syn> operator+(const cmatrix<T,Syn>& m0, const cmatrix<T,Syn>& m1) {
   return (cmatrix<T,Syn>(m0) += m1);
}

template <typename T, typename Syn>
cmatrix<T,Syn> operator-(const cmatrix<T,Syn>& m0, const cmatrix<T,Syn>& m1) {
   return (cmatrix<T,Syn>(m0) -= m1);
}

template <typename T, typename Syn>
cmatrix<T,Syn> operator*(const cmatrix<T,Syn>& m0, const cmatrix<T,Syn>& m1) {
   return (cmatrix<T,Syn>(m0) *= m1);
}

template <typename T, typename Syn>
cmatrix<T,Syn> operator/(const cmatrix<T,Syn>& m0, const cmatrix<T,Syn>& m1) {
   return (cmatrix<T,Syn>(m0) /= m1);
}

/***************************************************************************
 * Binary cmatrix-real comparators.
 ***************************************************************************/

template <typename T, typename Syn>
matrix<bool,Syn> operator==(const cmatrix<T,Syn>& m, const T& t) {
   return (m == complex<T>(t));
}

template <typename T, typename Syn>
matrix<bool,Syn> operator!=(const cmatrix<T,Syn>& m, const T& t) {
   return (m != complex<T>(t));
}

template <typename T, typename Syn>
matrix<bool,Syn> operator<(const cmatrix<T,Syn>& m, const T& t) {
   return (m < complex<T>(t));
}

template <typename T, typename Syn>
matrix<bool,Syn> operator>(const cmatrix<T,Syn>& m, const T& t) {
   return (m > complex<T>(t));
}

template <typename T, typename Syn>
matrix<bool,Syn> operator<=(const cmatrix<T,Syn>& m, const T& t) {
   return (m <= complex<T>(t));
}

template <typename T, typename Syn>
matrix<bool,Syn> operator>=(const cmatrix<T,Syn>& m, const T& t) {
   return (m >= complex<T>(t));
}

template <typename T, typename Syn>
matrix<bool,Syn> operator==(const T& t, const cmatrix<T,Syn>& m) {
   return (complex<T>(t) == m);
}

template <typename T, typename Syn>
matrix<bool,Syn> operator!=(const T& t, const cmatrix<T,Syn>& m) {
   return (complex<T>(t) != m);
}

template <typename T, typename Syn>
matrix<bool,Syn> operator<(const T& t, const cmatrix<T,Syn>& m) {
   return (complex<T>(t) < m);
}

template <typename T, typename Syn>
matrix<bool,Syn> operator>(const T& t, const cmatrix<T,Syn>& m) {
   return (complex<T>(t) > m);
}

template <typename T, typename Syn>
matrix<bool,Syn> operator<=(const T& t, const cmatrix<T,Syn>& m) {
   return (complex<T>(t) <= m);
}

template <typename T, typename Syn>
matrix<bool,Syn> operator>=(const T& t, const cmatrix<T,Syn>& m) {
   return (complex<T>(t) >= m);
}

/***************************************************************************
 * Binary cmatrix-complex comparators.
 ***************************************************************************/

template <typename T, typename Syn>
matrix<bool,Syn> operator==(const cmatrix<T,Syn>& m, const complex<T>& z) {
   return (static_cast<const matrix<complex<T>,Syn>&>(m) == z);
}

template <typename T, typename Syn>
matrix<bool,Syn> operator!=(const cmatrix<T,Syn>& m, const complex<T>& z) {
   return (static_cast<const matrix<complex<T>,Syn>&>(m) != z);
}

template <typename T, typename Syn>
matrix<bool,Syn> operator<(const cmatrix<T,Syn>& m, const complex<T>& z) {
   return (static_cast<const matrix<complex<T>,Syn>&>(m) < z);
}

template <typename T, typename Syn>
matrix<bool,Syn> operator>(const cmatrix<T,Syn>& m, const complex<T>& z) {
   return (static_cast<const matrix<complex<T>,Syn>&>(m) > z);
}

template <typename T, typename Syn>
matrix<bool,Syn> operator<=(const cmatrix<T,Syn>& m, const complex<T>& z) {
   return (static_cast<const matrix<complex<T>,Syn>&>(m) <= z);
}

template <typename T, typename Syn>
matrix<bool,Syn> operator>=(const cmatrix<T,Syn>& m, const complex<T>& z) {
   return (static_cast<const matrix<complex<T>,Syn>&>(m) >= z);
}

template <typename T, typename Syn>
matrix<bool,Syn> operator==(const complex<T>& z, const cmatrix<T,Syn>& m) {
   return (z == static_cast<const matrix<complex<T>,Syn>&>(m));
}

template <typename T, typename Syn>
matrix<bool,Syn> operator!=(const complex<T>& z, const cmatrix<T,Syn>& m) {
   return (z != static_cast<const matrix<complex<T>,Syn>&>(m));
}

template <typename T, typename Syn>
matrix<bool,Syn> operator<(const complex<T>& z, const cmatrix<T,Syn>& m) {
   return (z < static_cast<const matrix<complex<T>,Syn>&>(m));
}

template <typename T, typename Syn>
matrix<bool,Syn> operator>(const complex<T>& z, const cmatrix<T,Syn>& m) {
   return (z > static_cast<const matrix<complex<T>,Syn>&>(m));
}

template <typename T, typename Syn>
matrix<bool,Syn> operator<=(const complex<T>& z, const cmatrix<T,Syn>& m) {
   return (z <= static_cast<const matrix<complex<T>,Syn>&>(m));
}

template <typename T, typename Syn>
matrix<bool,Syn> operator>=(const complex<T>& z, const cmatrix<T,Syn>& m) {
   return (z >= static_cast<const matrix<complex<T>,Syn>&>(m));
}

/***************************************************************************
 * Binary cmatrix-matrix comparators.
 ***************************************************************************/

template <typename T, typename Syn>
matrix<bool,Syn> operator==(const cmatrix<T,Syn>& m0, const matrix<T,Syn>& m1) {
   auto_read_read_lock<const Syn> rrlock(m0, m1);
   matrix<T,Syn>::assert_dims_equal(m0._dims, m1._dims);
   matrix<bool,Syn> m(m0._dims);
   for (unsigned long n = 0; n < m0._size; n++)
      m[n] = (m0._data[n] == m1._data[n]);
   return m;
}

template <typename T, typename Syn>
matrix<bool,Syn> operator!=(const cmatrix<T,Syn>& m0, const matrix<T,Syn>& m1) {
   auto_read_read_lock<const Syn> rrlock(m0, m1);
   matrix<T,Syn>::assert_dims_equal(m0._dims, m1._dims);
   matrix<bool,Syn> m(m0._dims);
   for (unsigned long n = 0; n < m0._size; n++)
      m[n] = (m0._data[n] != m1._data[n]);
   return m;
}

template <typename T, typename Syn>
matrix<bool,Syn> operator<(const cmatrix<T,Syn>& m0, const matrix<T,Syn>& m1) {
   auto_read_read_lock<const Syn> rrlock(m0, m1);
   matrix<T,Syn>::assert_dims_equal(m0._dims, m1._dims);
   matrix<bool,Syn> m(m0._dims);
   for (unsigned long n = 0; n < m0._size; n++)
      m[n] = (m0._data[n] < m1._data[n]);
   return m;
}

template <typename T, typename Syn>
matrix<bool,Syn> operator>(const cmatrix<T,Syn>& m0, const matrix<T,Syn>& m1) {
   auto_read_read_lock<const Syn> rrlock(m0, m1);
   matrix<T,Syn>::assert_dims_equal(m0._dims, m1._dims);
   matrix<bool,Syn> m(m0._dims);
   for (unsigned long n = 0; n < m0._size; n++)
      m[n] = (m0._data[n] > m1._data[n]);
   return m;
}

template <typename T, typename Syn>
matrix<bool,Syn> operator<=(const cmatrix<T,Syn>& m0, const matrix<T,Syn>& m1) {
   auto_read_read_lock<const Syn> rrlock(m0, m1);
   matrix<T,Syn>::assert_dims_equal(m0._dims, m1._dims);
   matrix<bool,Syn> m(m0._dims);
   for (unsigned long n = 0; n < m0._size; n++)
      m[n] = (m0._data[n] <= m1._data[n]);
   return m;
}

template <typename T, typename Syn>
matrix<bool,Syn> operator>=(const cmatrix<T,Syn>& m0, const matrix<T,Syn>& m1) {
   auto_read_read_lock<const Syn> rrlock(m0, m1);
   matrix<T,Syn>::assert_dims_equal(m0._dims, m1._dims);
   matrix<bool,Syn> m(m0._dims);
   for (unsigned long n = 0; n < m0._size; n++)
      m[n] = (m0._data[n] >= m1._data[n]);
   return m;
}

template <typename T, typename Syn>
matrix<bool,Syn> operator==(const matrix<T,Syn>& m0, const cmatrix<T,Syn>& m1) {
   return (m1 == m0);
}

template <typename T, typename Syn>
matrix<bool,Syn> operator!=(const matrix<T,Syn>& m0, const cmatrix<T,Syn>& m1) {
   return (m1 != m0);
}

template <typename T, typename Syn>
matrix<bool,Syn> operator<(const matrix<T,Syn>& m0, const cmatrix<T,Syn>& m1) {
   return (m1 > m0);
}

template <typename T, typename Syn>
matrix<bool,Syn> operator>(const matrix<T,Syn>& m0, const cmatrix<T,Syn>& m1) {
   return (m1 < m0);
}

template <typename T, typename Syn>
matrix<bool,Syn> operator<=(const matrix<T,Syn>& m0, const cmatrix<T,Syn>& m1) {
   return (m1 >= m0);
}

template <typename T, typename Syn>
matrix<bool,Syn> operator>=(const matrix<T,Syn>& m0, const cmatrix<T,Syn>& m1) {
   return (m1 <= m0);
}

/***************************************************************************
 * Binary cmatrix-cmatrix comparators.
 ***************************************************************************/

template <typename T, typename Syn>
matrix<bool,Syn> operator==(const cmatrix<T,Syn>& m0, const cmatrix<T,Syn>& m1)
{
   return (
      static_cast<const matrix<complex<T>,Syn>&>(m0) == 
      static_cast<const matrix<complex<T>,Syn>&>(m1)
   );
}

template <typename T, typename Syn>
matrix<bool,Syn> operator!=(const cmatrix<T,Syn>& m0, const cmatrix<T,Syn>& m1)
{
   return (
      static_cast<const matrix<complex<T>,Syn>&>(m0) != 
      static_cast<const matrix<complex<T>,Syn>&>(m1)
   );
}

template <typename T, typename Syn>
matrix<bool,Syn> operator<(const cmatrix<T,Syn>& m0, const cmatrix<T,Syn>& m1)
{
   return (
      static_cast<const matrix<complex<T>,Syn>&>(m0) <
      static_cast<const matrix<complex<T>,Syn>&>(m1)
   );
}

template <typename T, typename Syn>
matrix<bool,Syn> operator>(const cmatrix<T,Syn>& m0, const cmatrix<T,Syn>& m1)
{
   return (
      static_cast<const matrix<complex<T>,Syn>&>(m0) >
      static_cast<const matrix<complex<T>,Syn>&>(m1)
   );
}

template <typename T, typename Syn>
matrix<bool,Syn> operator<=(const cmatrix<T,Syn>& m0, const cmatrix<T,Syn>& m1)
{
   return (
      static_cast<const matrix<complex<T>,Syn>&>(m0) <= 
      static_cast<const matrix<complex<T>,Syn>&>(m1)
   );
}

template <typename T, typename Syn>
matrix<bool,Syn> operator>=(const cmatrix<T,Syn>& m0, const cmatrix<T,Syn>& m1)
{
   return (
      static_cast<const matrix<complex<T>,Syn>&>(m0) >= 
      static_cast<const matrix<complex<T>,Syn>&>(m1)
   );
}

/***************************************************************************
 * Unary operators.
 ***************************************************************************/

/*
 * Unary arithmetic operators.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> cmatrix<T,Syn>::operator+() const {
   matrix<complex<T>,Syn> m = this->matrix<complex<T>,Syn>::operator+();
   return cmatrix<T,Syn>::to_cmatrix(m);
}

template <typename T, typename Syn>
cmatrix<T,Syn> cmatrix<T,Syn>::operator-() const {
   matrix<complex<T>,Syn> m = this->matrix<complex<T>,Syn>::operator-();
   return cmatrix<T,Syn>::to_cmatrix(m);
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::operator++() {
   this->matrix<complex<T>,Syn>::operator++();
   return *this;
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::operator--() {
   this->matrix<complex<T>,Syn>::operator--();
   return *this;
}

template <typename T, typename Syn>
cmatrix<T,Syn> cmatrix<T,Syn>::operator++(int) {
   matrix<complex<T>,Syn> m = this->matrix<complex<T>,Syn>::operator++(0);
   return cmatrix<T,Syn>::to_cmatrix(m);
}

template <typename T, typename Syn>
cmatrix<T,Syn> cmatrix<T,Syn>::operator--(int) {
   matrix<complex<T>,Syn> m = this->matrix<complex<T>,Syn>::operator--(0);
   return cmatrix<T,Syn>::to_cmatrix(m);
}

/*
 * Unary element multiplicative inversion operator.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> cmatrix<T,Syn>::operator~() const {
   matrix<complex<T>,Syn> m = this->matrix<complex<T>,Syn>::operator~();
   return cmatrix<T,Syn>::to_cmatrix(m);
}

/***************************************************************************
 * Transcendentals and other functions on matrix elements.
 ***************************************************************************/

template <typename T, typename Syn>
cmatrix<T,Syn> acos(const matrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   cmatrix<T,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result._data[n] = acos(complex<T>(m._data[n]));
   return m_result;
}

template <typename T, typename Syn>
cmatrix<T,Syn> asin(const matrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   cmatrix<T,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result._data[n] = asin(complex<T>(m._data[n]));
   return m_result;
}

template <typename T, typename Syn>
cmatrix<T,Syn> acosh(const matrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   cmatrix<T,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result._data[n] = acosh(complex<T>(m._data[n]));
   return m_result;
}

template <typename T, typename Syn>
cmatrix<T,Syn> atanh(const matrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   cmatrix<T,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result._data[n] = atanh(complex<T>(m._data[n]));
   return m_result;
}

template <typename T, typename Syn>
cmatrix<T,Syn> sqrt(const matrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   cmatrix<T,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result._data[n] = sqrt(complex<T>(m._data[n]));
   return m_result;
}

template <typename T, typename Syn>
cmatrix<T,Syn> log(const matrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   cmatrix<T,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result._data[n] = log(complex<T>(m._data[n]));
   return m_result;
}

template <typename T, typename Syn>
cmatrix<T,Syn> pow(const matrix<T,Syn>& m, const T& t) {
   const T t_copy(t);
   auto_read_lock<const Syn> rlock(m);
   cmatrix<T,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result._data[n] = pow(complex<T>(m._data[n]), t_copy);
   return m_result;
}

template <typename T, typename Syn>
cmatrix<T,Syn> pow(const matrix<T,Syn>& m, const complex<T>& z) {
   const complex<T> z_copy(z);
   auto_read_lock<const Syn> rlock(m);
   cmatrix<T,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result._data[n] = pow(m._data[n], z_copy);
   return m_result;
}

/***************************************************************************
 * Transcendentals and other functions on cmatrix elements.
 ***************************************************************************/

template <typename T, typename Syn>
cmatrix<T,Syn> cos(const cmatrix<T,Syn>& m) {
   matrix<complex<T>,Syn> m_result = cos(
      static_cast<const matrix<complex<T>,Syn>&>(m)
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn> sin(const cmatrix<T,Syn>& m) {
   matrix<complex<T>,Syn> m_result = sin(
      static_cast<const matrix<complex<T>,Syn>&>(m)
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn> tan(const cmatrix<T,Syn>& m) {
   matrix<complex<T>,Syn> m_result = tan(
      static_cast<const matrix<complex<T>,Syn>&>(m)
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn> cosh(const cmatrix<T,Syn>& m) {
   matrix<complex<T>,Syn> m_result = cosh(
      static_cast<const matrix<complex<T>,Syn>&>(m)
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn> sinh(const cmatrix<T,Syn>& m) {
   matrix<complex<T>,Syn> m_result = sinh(
      static_cast<const matrix<complex<T>,Syn>&>(m)
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn> tanh(const cmatrix<T,Syn>& m) {
   matrix<complex<T>,Syn> m_result = tanh(
      static_cast<const matrix<complex<T>,Syn>&>(m)
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn> acos(const cmatrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   cmatrix<T,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result._data[n] = acos(m._data[n]);
   return m_result;
}

template <typename T, typename Syn>
cmatrix<T,Syn> asin(const cmatrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   cmatrix<T,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result._data[n] = asin(m._data[n]);
   return m_result;
}

template <typename T, typename Syn>
cmatrix<T,Syn> atan(const cmatrix<T,Syn>& m) {
   matrix<complex<T>,Syn> m_result = atan(
      static_cast<const matrix<complex<T>,Syn>&>(m)
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn> acosh(const cmatrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   cmatrix<T,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result._data[n] = acosh(m._data[n]);
   return m_result;
}

template <typename T, typename Syn>
cmatrix<T,Syn> asinh(const cmatrix<T,Syn>& m) {
   matrix<complex<T>,Syn> m_result = asinh(
      static_cast<const matrix<complex<T>,Syn>&>(m)
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn> atanh(const cmatrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   cmatrix<T,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result._data[n] = atanh(m._data[n]);
   return m_result;
}

template <typename T, typename Syn>
cmatrix<T,Syn> sqrt(const cmatrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   cmatrix<T,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result._data[n] = sqrt(m._data[n]);
   return m_result;
}

template <typename T, typename Syn>
cmatrix<T,Syn> log(const cmatrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   cmatrix<T,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result._data[n] = log(m._data[n]);
   return m_result;
}

template <typename T, typename Syn>
cmatrix<T,Syn> exp(const cmatrix<T,Syn>& m) {
   matrix<complex<T>,Syn> m_result = exp(
      static_cast<const matrix<complex<T>,Syn>&>(m)
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn> pow(const cmatrix<T,Syn>& m, int p) {
   matrix<complex<T>,Syn> m_result = pow(
      static_cast<const matrix<complex<T>,Syn>&>(m), p
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn> pow(const cmatrix<T,Syn>& m, const T& t) {
   const T t_copy(t);
   auto_read_lock<const Syn> rlock(m);
   cmatrix<T,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result._data[n] = pow(m._data[n], t_copy);
   return m_result;
}

template <typename T, typename Syn>
cmatrix<T,Syn> pow(const cmatrix<T,Syn>& m, const complex<T>& z) {
   const complex<T> z_copy(z);
   auto_read_lock<const Syn> rlock(m);
   cmatrix<T,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result._data[n] = pow(m._data[n], z_copy);
   return m_result;
}

/*
 * Two-parameter arctangent of the real parts of the corresponding elements.
 */
template <typename T, typename Syn>
matrix<T,Syn> atan2(const cmatrix<T,Syn>& my, const cmatrix<T,Syn>& mx) {
   auto_read_read_lock<const Syn> rrlock(my, mx);
   matrix<T,Syn>::assert_dims_equal(my._dims, mx._dims);
   matrix<T,Syn> m_result(my._dims);
   for (unsigned long n = 0; n < my._size; n++)
      m_result._data[n] = atan2(my._data[n], mx._data[n]);
   return m_result;
}

/*
 * Magnitude of elements.
 */
template <typename T, typename Syn>
matrix<T,Syn> abs(const cmatrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   matrix<T,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result._data[n] = abs(m._data[n]);
   return m_result;
}

/*
 * Argument of elements.
 */
template <typename T, typename Syn>
matrix<T,Syn> arg(const cmatrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   matrix<T,Syn> m_result(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_result._data[n] = arg(m._data[n]);
   return m_result;
}

/***************************************************************************
 * Real/imaginary components.
 ***************************************************************************/

/*
 * Check if matrix is real-valued.
 */
template <typename T, typename Syn>
bool cmatrix<T,Syn>::is_real() const {
   auto_read_lock<const Syn> rlock(*this);
   bool flag = true;
   for (unsigned long n = 0; ((n < this->size) && (flag)); n++)
      flag = this->_data[n].is_real();
   return flag;
}

/*
 * Real component.
 */
template <typename T, typename Syn>
matrix<T,Syn> cmatrix<T,Syn>::real() const {
   auto_read_lock<const Syn> rlock(*this);
   matrix<T,Syn> m(this->_dims);
   for (unsigned long n = 0; n < this->_size; n++)
      m._data[n] = this->_data[n].real();
   return m;
}

/*
 * Imaginary component. 
 */
template <typename T, typename Syn>
matrix<T,Syn> cmatrix<T,Syn>::imag() const {
   auto_read_lock<const Syn> rlock(*this);
   matrix<T,Syn> m(this->_dims);
   for (unsigned long n = 0; n < this->_size; n++)
      m._data[n] = this->_data[n].imag();
   return m;
}

/*
 * Real and imaginary components.
 */
template <typename T, typename Syn>
void cmatrix<T,Syn>::parts(matrix<T,Syn>& m_real, matrix<T,Syn>& m_imag) const {
   /* copy real and imaginary components */
   matrix<T,Syn> m_real_copy;
   matrix<T,Syn> m_imag_copy;
   {
      auto_read_lock<const Syn> rlock(*this);
      m_real_copy.resize(this->_dims);
      m_imag_copy.resize(this->_dims);
      for (unsigned long n = 0; n < this->_size; n++)
         this->_data[n].parts(m_real_copy._data[n], m_imag_copy._data[n]);
   }
   /* assign real and imaginary components */
   auto_write_write_lock<const Syn> wwlock(m_real, m_imag);
   matrix<T,Syn>::swap(m_real, m_real_copy);
   matrix<T,Syn>::swap(m_imag, m_imag_copy);
}

/*
 * Complex conjugate.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> conj(const cmatrix<T,Syn>& m) {
   auto_read_lock<const Syn> rlock(m);
   cmatrix<T,Syn> m_conj(m._dims);
   for (unsigned long n = 0; n < m._size; n++)
      m_conj._data[n] = conj(m._data[n]);
   return m_conj;
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::conj() {
   auto_write_lock<const Syn> wlock(*this);
   for (unsigned long n = 0; n < this->_size; n++)
      this->_data[n] = conj(this->_data[n]);
   return *this;
}

/***************************************************************************
 * Convolution, dot product, and element product.
 ***************************************************************************/

/*
 * Convolution.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> conv(
   const cmatrix<T,Syn>& m0, const cmatrix<T,Syn>& m1)
{
   matrix<complex<T>,Syn> m = conv(
      static_cast<const matrix<complex<T>,Syn>&>(m0),
      static_cast<const matrix<complex<T>,Syn>&>(m1)
   );
   return cmatrix<T,Syn>::to_cmatrix(m);
}

/*
 * Convolution.
 * Crop the result to be no larger than the left input matrix.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> conv_crop(
   const cmatrix<T,Syn>& m0, const cmatrix<T,Syn>& m1)
{
   matrix<complex<T>,Syn> m = conv_crop(
      static_cast<const matrix<complex<T>,Syn>&>(m0),
      static_cast<const matrix<complex<T>,Syn>&>(m1)
   );
   return cmatrix<T,Syn>::to_cmatrix(m);
}

/*
 * Convolution.
 * Crop the result to be no larger than the left input matrix and to
 * include only the central portion for which no zero padding was required.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> conv_crop_strict(
   const cmatrix<T,Syn>& m0, const cmatrix<T,Syn>& m1)
{
   matrix<complex<T>,Syn> m = conv_crop_strict(
      static_cast<const matrix<complex<T>,Syn>&>(m0),
      static_cast<const matrix<complex<T>,Syn>&>(m1)
   );
   return cmatrix<T,Syn>::to_cmatrix(m);
}

/*
 * Dot product.
 */
template <typename T, typename Syn>
complex<T> dot(const cmatrix<T,Syn>& m0, const cmatrix<T,Syn>& m1) {
   return dot(
      static_cast<const matrix<complex<T>,Syn>&>(m0),
      static_cast<const matrix<complex<T>,Syn>&>(m1)
   );
}

/*
 * Element product.
 * Compute the product of corresponding matrix elements.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> prod(
   const cmatrix<T,Syn>& m0, const cmatrix<T,Syn>& m1)
{
   matrix<complex<T>,Syn> m = prod(
      static_cast<const matrix<complex<T>,Syn>&>(m0),
      static_cast<const matrix<complex<T>,Syn>&>(m1)
   );
   return cmatrix<T,Syn>::to_cmatrix(m);
}

/***************************************************************************
 * Minimum, maximum, sum, product.
 ***************************************************************************/

/*
 * Minimum, maximum, sum, product of all matrix elements.
 * Note that min/max operate on lexicographic ordering.
 */
template <typename T, typename Syn>
complex<T> min(const cmatrix<T,Syn>& m) {
   return min(static_cast<const matrix<complex<T>,Syn>&>(m));
}

template <typename T, typename Syn>
complex<T> max(const cmatrix<T,Syn>& m) {
   return max(static_cast<const matrix<complex<T>,Syn>&>(m));
}

template <typename T, typename Syn>
complex<T> sum(const cmatrix<T,Syn>& m) {
   return sum(static_cast<const matrix<complex<T>,Syn>&>(m));
}

template <typename T, typename Syn>
complex<T> prod(const cmatrix<T,Syn>& m) {
   return prod(static_cast<const matrix<complex<T>,Syn>&>(m));
}
 
/*
 * Minimum, maximum, sum, product along specified dimension.
 * Note that min/max operate on lexicographic ordering.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> min(const cmatrix<T,Syn>& m, unsigned long d) {
   matrix<complex<T>,Syn> m_result = min(
      static_cast<const matrix<complex<T>,Syn>&>(m), d
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn> max(const cmatrix<T,Syn>& m, unsigned long d) {
   matrix<complex<T>,Syn> m_result = max(
      static_cast<const matrix<complex<T>,Syn>&>(m), d
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn> sum(const cmatrix<T,Syn>& m, unsigned long d) {
   matrix<complex<T>,Syn> m_result = sum(
      static_cast<const matrix<complex<T>,Syn>&>(m), d
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn> prod(const cmatrix<T,Syn>& m, unsigned long d) {
   matrix<complex<T>,Syn> m_result = prod(
      static_cast<const matrix<complex<T>,Syn>&>(m), d
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

/***************************************************************************
 * Cumulative sum, product.
 ***************************************************************************/

/*
 * Cumulative sum, product of all matrix elements.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> cumsum(const cmatrix<T,Syn>& m) {
   matrix<complex<T>,Syn> m_result = cumsum(
      static_cast<const matrix<complex<T>,Syn>&>(m)
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn> cumprod(const cmatrix<T,Syn>& m) {
   matrix<complex<T>,Syn> m_result = cumprod(
      static_cast<const matrix<complex<T>,Syn>&>(m)
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

/*
 * Cumulative sum, product along specified dimension.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> cumsum(const cmatrix<T,Syn>& m, unsigned long d) {
   matrix<complex<T>,Syn> m_result = cumsum(
      static_cast<const matrix<complex<T>,Syn>&>(m), d
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn> cumprod(const cmatrix<T,Syn>& m, unsigned long d) {
   matrix<complex<T>,Syn> m_result = cumprod(
      static_cast<const matrix<complex<T>,Syn>&>(m), d
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

/***************************************************************************
 * Mean, variance.
 ***************************************************************************/

/*
 * Mean, variance of all matrix elements.
 */
template <typename T, typename Syn>
complex<T> mean(const cmatrix<T,Syn>& m) {
   return mean(static_cast<const matrix<complex<T>,Syn>&>(m));
}

template <typename T, typename Syn>
complex<T> var(const cmatrix<T,Syn>& m) {
   return var(static_cast<const matrix<complex<T>,Syn>&>(m));
}

/*
 * Mean, variance along specified dimension.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> mean(const cmatrix<T,Syn>& m, unsigned long d) {
   matrix<complex<T>,Syn> m_result = mean(
      static_cast<const matrix<complex<T>,Syn>&>(m), d
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn> var(const cmatrix<T,Syn>& m, unsigned long d) {
   matrix<complex<T>,Syn> m_result = var(
      static_cast<const matrix<complex<T>,Syn>&>(m), d
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
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
cmatrix<T,Syn> gradient(
   const cmatrix<T,Syn>& m, unsigned long d)
{
   matrix<complex<T>,Syn> m_result = gradient(
      static_cast<const matrix<complex<T>,Syn>&>(m), d
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}


template <typename T, typename Syn>
cmatrix<T,Syn> gradient(
   const cmatrix<T,Syn>& m, unsigned long d, const complex<T>& z)
{
   matrix<complex<T>,Syn> m_result = gradient(
      static_cast<const matrix<complex<T>,Syn>&>(m), d, z
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

/***************************************************************************
 * Find.
 ***************************************************************************/

/*
 * Find locations of elements with the specified value.
 */
template <typename T, typename Syn>
array<unsigned long> cmatrix<T,Syn>::find_value(const T& t) const {
   return this->find_value(complex<T>(t));
}

/***************************************************************************
 * Diagonal.
 ***************************************************************************/

/*
 * Diagonal entries (returns diagonal of matrix in a vector).
 */
template <typename T, typename Syn>
cmatrix<T,Syn> diag(const cmatrix<T,Syn>& m) {
   return m.diag();
}

template <typename T, typename Syn>
cmatrix<T,Syn> cmatrix<T,Syn>::diag() const {
   matrix<complex<T>,Syn> dgnl = this->matrix<complex<T>,Syn>::diag();
   return cmatrix<T,Syn>::to_cmatrix(dgnl);
}

/***************************************************************************
 * Inverse and row reduce.
 ***************************************************************************/

/*
 * Inverse.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> inv(const cmatrix<T,Syn>& m, const T& tol) {
   matrix<complex<T>,Syn> m_result = inv(
      static_cast<const matrix<complex<T>,Syn>&>(m), complex<T>(tol)
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::inv(const T& tol) {
   this->matrix<complex<T>,Syn>::inv(complex<T>(tol));
   return *this;
}

/*
 * Row echelon form.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> ref(const cmatrix<T,Syn>& m, const T& tol) {
   matrix<complex<T>,Syn> m_result = ref(
      static_cast<const matrix<complex<T>,Syn>&>(m), complex<T>(tol)
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::ref(const T& tol) {
   this->matrix<complex<T>,Syn>::ref(complex<T>(tol));
   return *this;
}

/*
 * Reduced row echelon form.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> rref(const cmatrix<T,Syn>& m, const T& tol) {
   matrix<complex<T>,Syn> m_result = rref(
      static_cast<const matrix<complex<T>,Syn>&>(m), complex<T>(tol)
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::rref(const T& tol) {
   this->matrix<complex<T>,Syn>::rref(complex<T>(tol));
   return *this;
}

/***************************************************************************
 * Determinant and rank.
 ***************************************************************************/

/*
 * Determinant.
 */
template <typename T, typename Syn>
complex<T> cmatrix<T,Syn>::det(const T& tol) const {
   return this->matrix<complex<T>,Syn>::det(complex<T>(tol));
}

/*
 * Rank.
 */
template <typename T, typename Syn>
unsigned long cmatrix<T,Syn>::rank(const T& tol) const {
   return this->matrix<complex<T>,Syn>::rank(complex<T>(tol));
}

/***************************************************************************
 * Dimensions and size.
 ***************************************************************************/

/* 
 * Transpose.
 * Reverse dimension order.
 * The input matrix is regarded as being at least two-dimensional.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> transpose(const cmatrix<T,Syn>& m) {
   matrix<complex<T>,Syn> m_result = transpose(
      static_cast<const matrix<complex<T>,Syn>&>(m)
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::transpose() {
   this->matrix<complex<T>,Syn>::transpose();
   return *this;
}

/*
 * Dimension shifting.
 * Shift dimensions of the matrix to the left, wrapping the leading
 * dimensions around to the right.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> shift_dimensions(
   const cmatrix<T,Syn>& m,
   unsigned long n_shift)
{
   matrix<complex<T>,Syn> m_result = shift_dimensions(
      static_cast<const matrix<complex<T>,Syn>&>(m), n_shift
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::shift_dimensions(
   unsigned long n_shift)
{
   this->matrix<complex<T>,Syn>::shift_dimensions(n_shift);
   return *this;
}

/*
 * Dimension permutation.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> permute_dimensions(
   const cmatrix<T,Syn>& m,
   const array<unsigned long>& order)
{
   matrix<complex<T>,Syn> m_result = permute_dimensions(
      static_cast<const matrix<complex<T>,Syn>&>(m), order
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::permute_dimensions(
   const array<unsigned long>& order)
{
   this->matrix<complex<T>,Syn>::permute_dimensions(order);
   return *this;
}

/*
 * Singleton dimension elimination.
 * Remove singleton dimensions, leaving the matrix elements unchanged.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> squeeze_dimensions(const cmatrix<T,Syn>& m) {
   matrix<complex<T>,Syn> m_result = squeeze_dimensions(
      static_cast<const matrix<complex<T>,Syn>&>(m)
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::squeeze_dimensions() {
   this->matrix<complex<T>,Syn>::squeeze_dimensions();
   return *this;
}

/***************************************************************************
 * Resize.
 ***************************************************************************/

/*
 * Resize (vector).
 * Return the resized vector.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> resize(
   const cmatrix<T,Syn>& m, unsigned long length)
{
   matrix<complex<T>,Syn> m_result = resize(
      static_cast<const matrix<complex<T>,Syn>&>(m), length
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn> resize(
   const cmatrix<T,Syn>& m, unsigned long length, const T& t)
{
   matrix<complex<T>,Syn> m_result = resize(
      static_cast<const matrix<complex<T>,Syn>&>(m), length, complex<T>(t)
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn> resize(
   const cmatrix<T,Syn>& m, unsigned long length, const complex<T>& z)
{
   matrix<complex<T>,Syn> m_result = resize(
      static_cast<const matrix<complex<T>,Syn>&>(m), length, z
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::resize(
   unsigned long length)
{
   this->matrix<complex<T>,Syn>::resize(length);
   return *this;
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::resize(
   unsigned long length, const T& t)
{
   this->matrix<complex<T>,Syn>::resize(length, complex<T>(t));
   return *this;
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::resize(
   unsigned long length, const complex<T>& z)
{
   this->matrix<complex<T>,Syn>::resize(length, z);
   return *this;
}

/*
 * Resize (2D matrix).
 * Return the resized M x N matrix.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> resize(
   const cmatrix<T,Syn>& m, unsigned long M, unsigned long N)
{
   matrix<complex<T>,Syn> m_result = resize(
      static_cast<const matrix<complex<T>,Syn>&>(m), M, N
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn> resize(
   const cmatrix<T,Syn>& m, unsigned long M, unsigned long N, const T& t)
{
   matrix<complex<T>,Syn> m_result = resize(
      static_cast<const matrix<complex<T>,Syn>&>(m), M, N, complex<T>(t)
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn> resize(
   const cmatrix<T,Syn>& m, unsigned long M, unsigned long N, const complex<T>& z)
{
   matrix<complex<T>,Syn> m_result = resize(
      static_cast<const matrix<complex<T>,Syn>&>(m), M, N, z
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::resize(
   unsigned long M, unsigned long N)
{
   this->matrix<complex<T>,Syn>::resize(M, N);
   return *this;
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::resize(
   unsigned long M, unsigned long N, const T& t)
{
   this->matrix<complex<T>,Syn>::resize(M, N, complex<T>(t));
   return *this;
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::resize(
   unsigned long M, unsigned long N, const complex<T>& z)
{
   this->matrix<complex<T>,Syn>::resize(M, N, z);
   return *this;
}

/*
 * Resize (multi-dimensional matrix).
 * Return the resized matrix.
 * The dimensionality of the matrix must not change.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> resize(
   const cmatrix<T,Syn>& m, const array<unsigned long>& dims)
{
   matrix<complex<T>,Syn> m_result = resize(
      static_cast<const matrix<complex<T>,Syn>&>(m), dims
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn> resize(
   const cmatrix<T,Syn>& m, const array<unsigned long>& dims, const T& t)
{
   matrix<complex<T>,Syn> m_result = resize(
      static_cast<const matrix<complex<T>,Syn>&>(m), dims, complex<T>(t)
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn> resize(
   const cmatrix<T,Syn>& m, const array<unsigned long>& dims, const complex<T>& z)
{
   matrix<complex<T>,Syn> m_result = resize(
      static_cast<const matrix<complex<T>,Syn>&>(m), dims, z
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::resize(
   const array<unsigned long>& dims)
{
   this->matrix<complex<T>,Syn>::resize(dims);
   return *this;
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::resize(
   const array<unsigned long>& dims, const T& t)
{
   this->matrix<complex<T>,Syn>::resize(dims, complex<T>(t));
   return *this;
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::resize(
   const array<unsigned long>& dims, const complex<T>& z)
{
   this->matrix<complex<T>,Syn>::resize(dims, z);
   return *this;
}

/***************************************************************************
 * Reshape.
 ***************************************************************************/

/*
 * Reshape into a vector.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> vector(const cmatrix<T,Syn>& m) {
   matrix<complex<T>,Syn> m_result = vector(
      static_cast<const matrix<complex<T>,Syn>&>(m)
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::vector() {
   this->matrix<complex<T>,Syn>::vector();
   return *this;
}

/*
 * Reshape into 2D M x N matrix.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> reshape(
   const cmatrix<T,Syn>& m, unsigned long M, unsigned long N)
{
   matrix<complex<T>,Syn> m_result = reshape(
      static_cast<const matrix<complex<T>,Syn>&>(m), M, N
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::reshape(unsigned long M, unsigned long N) {
   this->matrix<complex<T>,Syn>::reshape(M, N);
   return *this;
}

/*
 * Reshape into multi-dimensional matrix.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> reshape(
   const cmatrix<T,Syn>& m, const array<unsigned long>& dims)
{
   matrix<complex<T>,Syn> m_result = reshape(
      static_cast<const matrix<complex<T>,Syn>&>(m), dims
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::reshape(const array<unsigned long>& dims) {
   this->matrix<complex<T>,Syn>::reshape(dims);
   return *this;
}

/***************************************************************************
 * Replication.
 ***************************************************************************/

/*
 * Replicate matrix.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> repmat(
   const cmatrix<T,Syn>& m, unsigned long M, unsigned long N)
{
   matrix<complex<T>,Syn> m_result = repmat(
      static_cast<const matrix<complex<T>,Syn>&>(m), M, N
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}
   
template <typename T, typename Syn>
cmatrix<T,Syn> repmat(
   const cmatrix<T,Syn>& m, const array<unsigned long>& dims)
{
   matrix<complex<T>,Syn> m_result = repmat(
      static_cast<const matrix<complex<T>,Syn>&>(m), dims
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

/***************************************************************************
 * Concatenation.
 ***************************************************************************/

/*
 * Concatenate matrices.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> vertcat(
   const cmatrix<T,Syn>& m_top, const cmatrix<T,Syn>& m_bottom)
{
   matrix<complex<T>,Syn> m_result = vertcat(
      static_cast<const matrix<complex<T>,Syn>&>(m_top),
      static_cast<const matrix<complex<T>,Syn>&>(m_bottom)
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}
   
template <typename T, typename Syn>
cmatrix<T,Syn> horzcat(
   const cmatrix<T,Syn>& m_left, const cmatrix<T,Syn>& m_right)
{
   matrix<complex<T>,Syn> m_result = horzcat(
      static_cast<const matrix<complex<T>,Syn>&>(m_left),
      static_cast<const matrix<complex<T>,Syn>&>(m_right)
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

/*
 * Concatenate matrices along the specified dimension.
 */
template <typename T, typename Syn>
cmatrix<T,Syn> concat(
   const cmatrix<T,Syn>& m0, const cmatrix<T,Syn>& m1, unsigned long d)
{
   matrix<complex<T>,Syn> m_result = concat(
      static_cast<const matrix<complex<T>,Syn>&>(m0),
      static_cast<const matrix<complex<T>,Syn>&>(m1),
      d
   );
   return cmatrix<T,Syn>::to_cmatrix(m_result);
}

/***************************************************************************
 * Reverse.
 ***************************************************************************/

/*
 * Reverse element order (treating the matrix as a linear array).
 */
template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::reverse() {
   this->matrix<complex<T>,Syn>::reverse();
   return *this;
}

/*
 * Reverse element order along the specified dimension.
 */
template <typename T, typename Syn>
cmatrix<T,Syn>& cmatrix<T,Syn>::reverse(unsigned long d) {
   this->matrix<complex<T>,Syn>::reverse(d);
   return *this;
}

} /* namespace matrices */
} /* namespace math */

#endif
