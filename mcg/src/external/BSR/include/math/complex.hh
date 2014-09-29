/*
 * Complex numbers (thread-safe).
 */
#ifndef MATH__COMPLEX_HH
#define MATH__COMPLEX_HH

#include "concurrent/threads/synchronization/locks/auto_read_lock.hh"
#include "concurrent/threads/synchronization/locks/auto_read_read_lock.hh"
#include "concurrent/threads/synchronization/locks/auto_read_write_lock.hh"
#include "concurrent/threads/synchronization/locks/auto_write_lock.hh"
#include "concurrent/threads/synchronization/synchronizables/synchronizable.hh"
#include "concurrent/threads/synchronization/synchronizables/unsynchronized.hh"
#include "io/serialization/serial_input_stream.hh"
#include "io/serialization/serial_output_stream.hh"
#include "io/serialization/serializers.hh"
#include "io/streams/ostream.hh"
#include "lang/pointers/auto_ptr.hh"
#include "math/math.hh"

namespace math {
/*
 * Imports.
 */
using concurrent::threads::synchronization::locks::auto_read_lock;
using concurrent::threads::synchronization::locks::auto_read_read_lock;
using concurrent::threads::synchronization::locks::auto_read_write_lock;
using concurrent::threads::synchronization::locks::auto_write_lock;
using concurrent::threads::synchronization::synchronizables::synchronizable;
using concurrent::threads::synchronization::synchronizables::unsynchronized;
using io::serialization::serial_input_stream;
using io::serialization::serial_output_stream;
using io::serialization::serializer;
using io::serialization::serializers;
using io::streams::ostream;
using lang::pointers::auto_ptr;

/*
 * Declare existence of complex class.
 */
template <typename T, typename Syn> class complex;

/*
 * Declare prototypes for complex template friend functions.
 */
template <typename T, typename Syn>
ostream& operator<<(ostream&, const complex<T,Syn>&);

template <typename T, typename Syn>
complex<T,Syn> operator+(const complex<T,Syn>&, const T&);
template <typename T, typename Syn>
complex<T,Syn> operator-(const complex<T,Syn>&, const T&);
template <typename T, typename Syn>
complex<T,Syn> operator*(const complex<T,Syn>&, const T&);
template <typename T, typename Syn>
complex<T,Syn> operator/(const complex<T,Syn>&, const T&);

template <typename T, typename Syn>
complex<T,Syn> operator+(const T&, const complex<T,Syn>&);
template <typename T, typename Syn>
complex<T,Syn> operator-(const T&, const complex<T,Syn>&);
template <typename T, typename Syn>
complex<T,Syn> operator*(const T&, const complex<T,Syn>&);
template <typename T, typename Syn>
complex<T,Syn> operator/(const T&, const complex<T,Syn>&);

template <typename T, typename Syn>
complex<T,Syn> operator+(const complex<T,Syn>&, const complex<T,Syn>&);
template <typename T, typename Syn>
complex<T,Syn> operator-(const complex<T,Syn>&, const complex<T,Syn>&);
template <typename T, typename Syn>
complex<T,Syn> operator*(const complex<T,Syn>&, const complex<T,Syn>&);
template <typename T, typename Syn>
complex<T,Syn> operator/(const complex<T,Syn>&, const complex<T,Syn>&);

template <typename T, typename Syn>
bool operator==(const complex<T,Syn>&, const T&);
template <typename T, typename Syn>
bool operator!=(const complex<T,Syn>&, const T&);
template <typename T, typename Syn>
bool operator< (const complex<T,Syn>&, const T&);
template <typename T, typename Syn>
bool operator> (const complex<T,Syn>&, const T&);
template <typename T, typename Syn>
bool operator<=(const complex<T,Syn>&, const T&);
template <typename T, typename Syn>
bool operator>=(const complex<T,Syn>&, const T&);

template <typename T, typename Syn>
bool operator==(const T&, const complex<T,Syn>&);
template <typename T, typename Syn>
bool operator!=(const T&, const complex<T,Syn>&);
template <typename T, typename Syn>
bool operator< (const T&, const complex<T,Syn>&);
template <typename T, typename Syn>
bool operator> (const T&, const complex<T,Syn>&);
template <typename T, typename Syn>
bool operator<=(const T&, const complex<T,Syn>&);
template <typename T, typename Syn>
bool operator>=(const T&, const complex<T,Syn>&);

template <typename T, typename Syn>
bool operator==(const complex<T,Syn>&, const complex<T,Syn>&);
template <typename T, typename Syn>
bool operator!=(const complex<T,Syn>&, const complex<T,Syn>&);
template <typename T, typename Syn>
bool operator< (const complex<T,Syn>&, const complex<T,Syn>&);
template <typename T, typename Syn>
bool operator> (const complex<T,Syn>&, const complex<T,Syn>&);
template <typename T, typename Syn>
bool operator<=(const complex<T,Syn>&, const complex<T,Syn>&);
template <typename T, typename Syn>
bool operator>=(const complex<T,Syn>&, const complex<T,Syn>&);

template <typename T, typename Syn>
complex<T,Syn> cos(const complex<T,Syn>&);
template <typename T, typename Syn>
complex<T,Syn> sin(const complex<T,Syn>&);
template <typename T, typename Syn>
complex<T,Syn> tan(const complex<T,Syn>&);
template <typename T, typename Syn>
complex<T,Syn> cosh(const complex<T,Syn>&);
template <typename T, typename Syn>
complex<T,Syn> sinh(const complex<T,Syn>&);
template <typename T, typename Syn>
complex<T,Syn> tanh(const complex<T,Syn>&);

template <typename T, typename Syn>
complex<T,Syn> acos(const complex<T,Syn>&);
template <typename T, typename Syn>
complex<T,Syn> asin(const complex<T,Syn>&);
template <typename T, typename Syn>
complex<T,Syn> atan(const complex<T,Syn>&);
template <typename T, typename Syn>
complex<T,Syn> acosh(const complex<T,Syn>&);
template <typename T, typename Syn>
complex<T,Syn> asinh(const complex<T,Syn>&);
template <typename T, typename Syn>
complex<T,Syn> atanh(const complex<T,Syn>&);

template <typename T, typename Syn>
complex<T,Syn> sqrt(const complex<T,Syn>&);
template <typename T, typename Syn>
complex<T,Syn> log(const complex<T,Syn>&);
template <typename T, typename Syn>
complex<T,Syn> exp(const complex<T,Syn>&);
template <typename T, typename Syn>
complex<T,Syn> pow(const T&, const complex<T,Syn>&);
template <typename T, typename Syn>
complex<T,Syn> pow(const complex<T,Syn>&, int);
template <typename T, typename Syn>
complex<T,Syn> pow(const complex<T,Syn>&, const T&);
template <typename T, typename Syn>
complex<T,Syn> pow(const complex<T,Syn>&, const complex<T,Syn>&);

template <typename T, typename Syn>
complex<T,Syn> atan2(const complex<T,Syn>&, const complex<T,Syn>&);

template <typename T, typename Syn> T abs(const complex<T,Syn>&);
template <typename T, typename Syn> T arg(const complex<T,Syn>&);

template <typename T, typename Syn> complex<T,Syn> conj(const complex<T,Syn>&);

/*
 * Complex class.
 */
template <typename T = double, typename Syn = unsynchronized>
class complex : protected Syn {
public:
   /*
    * Friend classes.
    */
   template <typename U, typename S> friend class complex;

   /*
    * Constructors.
    */
   complex();
   explicit complex(const T&);
   explicit complex(const T&, const T&);

   /*
    * Copy constructors.
    */
   complex(const complex<T,Syn>&);
   template <typename U, typename S> complex(const complex<U,S>&);
    
   /*
    * Destructor.
    */
   virtual ~complex();
   
   /*
    * Polar form.
    */
   static complex<T,Syn> polar(const T& /* r */, const T& /* theta */);

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
   static auto_ptr< complex<T,Syn> > deserialize(
      serial_input_stream&,
      const serializer<T>& = serializers<T>::s_default()
   );

   /*
    * Formatted output to stream.
    */
   friend ostream& operator<< <T,Syn>(ostream&, const complex<T,Syn>&);

   /*
    * Assignment operators: complex-real.
    */
   complex<T,Syn>& operator=(const T&);
   complex<T,Syn>& operator+=(const T&);
   complex<T,Syn>& operator-=(const T&);
   complex<T,Syn>& operator*=(const T&);
   complex<T,Syn>& operator/=(const T&);
  
   /*
    * Assignment operators: complex-complex.
    */
   complex<T,Syn>& operator=(const complex<T,Syn>&);
   complex<T,Syn>& operator+=(const complex<T,Syn>&);
   complex<T,Syn>& operator-=(const complex<T,Syn>&);
   complex<T,Syn>& operator*=(const complex<T,Syn>&);
   complex<T,Syn>& operator/=(const complex<T,Syn>&);
   
   template <typename U, typename S> complex<T,Syn>&
      operator=(const complex<U,S>&);
   template <typename U, typename S> complex<T,Syn>&
      operator+=(const complex<U,S>&);
   template <typename U, typename S> complex<T,Syn>&
      operator-=(const complex<U,S>&);
   template <typename U, typename S> complex<T,Syn>&
      operator*=(const complex<U,S>&);
   template <typename U, typename S> complex<T,Syn>&
      operator/=(const complex<U,S>&);
   
   /*
    * Binary complex-real operators.
    */
   friend complex<T,Syn> operator+ <T,Syn>(const complex<T,Syn>&, const T&);
   friend complex<T,Syn> operator- <T,Syn>(const complex<T,Syn>&, const T&);
   friend complex<T,Syn> operator* <T,Syn>(const complex<T,Syn>&, const T&);
   friend complex<T,Syn> operator/ <T,Syn>(const complex<T,Syn>&, const T&);

   friend complex<T,Syn> operator+ <T,Syn>(const T&, const complex<T,Syn>&);
   friend complex<T,Syn> operator- <T,Syn>(const T&, const complex<T,Syn>&);
   friend complex<T,Syn> operator* <T,Syn>(const T&, const complex<T,Syn>&);
   friend complex<T,Syn> operator/ <T,Syn>(const T&, const complex<T,Syn>&);
   
   /*
    * Binary complex-complex operators.
    */
   friend complex<T,Syn>
      operator+ <T,Syn>(const complex<T,Syn>&, const complex<T,Syn>&);
   friend complex<T,Syn>
      operator- <T,Syn>(const complex<T,Syn>&, const complex<T,Syn>&);
   friend complex<T,Syn>
      operator* <T,Syn>(const complex<T,Syn>&, const complex<T,Syn>&);
   friend complex<T,Syn>
      operator/ <T,Syn>(const complex<T,Syn>&, const complex<T,Syn>&);

   /*
    * Binary complex-real comparators.
    */
   friend bool operator== <T,Syn>(const complex<T,Syn>&, const T&);
   friend bool operator!= <T,Syn>(const complex<T,Syn>&, const T&);
   friend bool operator<  <T,Syn>(const complex<T,Syn>&, const T&);
   friend bool operator>  <T,Syn>(const complex<T,Syn>&, const T&);
   friend bool operator<= <T,Syn>(const complex<T,Syn>&, const T&);
   friend bool operator>= <T,Syn>(const complex<T,Syn>&, const T&);
   
   friend bool operator== <T,Syn>(const T&, const complex<T,Syn>&);
   friend bool operator!= <T,Syn>(const T&, const complex<T,Syn>&);
   friend bool operator<  <T,Syn>(const T&, const complex<T,Syn>&);
   friend bool operator>  <T,Syn>(const T&, const complex<T,Syn>&);
   friend bool operator<= <T,Syn>(const T&, const complex<T,Syn>&);
   friend bool operator>= <T,Syn>(const T&, const complex<T,Syn>&);
   
   /*
    * Binary complex-complex comparators.
    */
   friend bool operator== <T,Syn>(const complex<T,Syn>&, const complex<T,Syn>&);
   friend bool operator!= <T,Syn>(const complex<T,Syn>&, const complex<T,Syn>&);
   friend bool operator<  <T,Syn>(const complex<T,Syn>&, const complex<T,Syn>&);
   friend bool operator>  <T,Syn>(const complex<T,Syn>&, const complex<T,Syn>&);
   friend bool operator<= <T,Syn>(const complex<T,Syn>&, const complex<T,Syn>&);
   friend bool operator>= <T,Syn>(const complex<T,Syn>&, const complex<T,Syn>&);
   
   /*
    * Unary arithmetic operators.
    */
   complex<T,Syn> operator+() const;
   complex<T,Syn> operator-() const;
   complex<T,Syn>& operator++();
   complex<T,Syn>& operator--();
   complex<T,Syn> operator++(int);
   complex<T,Syn> operator--(int);

   /*
    * Unary logical inversion operator.
    */
   inline bool operator!() const;
   
   /*
    * Unary multiplicative inversion operator.
    */
   complex<T,Syn> operator~() const;
  
   /*
    * Transcendentals.
    */
   friend complex<T,Syn> cos<T,Syn>(const complex<T,Syn>&);
   friend complex<T,Syn> sin<T,Syn>(const complex<T,Syn>&);
   friend complex<T,Syn> tan<T,Syn>(const complex<T,Syn>&);
   friend complex<T,Syn> cosh<T,Syn>(const complex<T,Syn>&);
   friend complex<T,Syn> sinh<T,Syn>(const complex<T,Syn>&);
   friend complex<T,Syn> tanh<T,Syn>(const complex<T,Syn>&);
   
   friend complex<T,Syn> acos<T,Syn>(const complex<T,Syn>&);
   friend complex<T,Syn> asin<T,Syn>(const complex<T,Syn>&);
   friend complex<T,Syn> atan<T,Syn>(const complex<T,Syn>&);
   friend complex<T,Syn> acosh<T,Syn>(const complex<T,Syn>&);
   friend complex<T,Syn> asinh<T,Syn>(const complex<T,Syn>&);
   friend complex<T,Syn> atanh<T,Syn>(const complex<T,Syn>&);
   
   friend complex<T,Syn> sqrt<T,Syn>(const complex<T,Syn>&);
   friend complex<T,Syn> log<T,Syn>(const complex<T,Syn>&);
   friend complex<T,Syn> exp<T,Syn>(const complex<T,Syn>&);
   friend complex<T,Syn> pow<T,Syn>(const T&, const complex<T,Syn>&);
   friend complex<T,Syn> pow<T,Syn>(const complex<T,Syn>&, int);
   friend complex<T,Syn> pow<T,Syn>(const complex<T,Syn>&, const T&);
   friend complex<T,Syn> pow<T,Syn>(
      const complex<T,Syn>&, const complex<T,Syn>&
   );
  
   /*
    * Two-parameter arctangent of the real parts.
    */
   friend complex<T,Syn> atan2<T,Syn>(
      const complex<T,Syn>& /* y */, const complex<T,Syn>& /* x */
   );

   /*
    * Magnitude.
    */
   friend T abs<T,Syn>(const complex<T,Syn>&);

   /*
    * Argument.
    */
   friend T arg<T,Syn>(const complex<T,Syn>&);

   /*
    * Check if real-valued (zero imaginary component).
    */
   bool is_real() const;

   /*
    * Check if imaginary-valued (zero real component).
    */
   bool is_imag() const;

   /*
    * Real component.
    */
   T real() const;
   
   /*
    * Imaginary component. 
    */
   T imag() const;
   
   /*
    * Real and imaginary components.
    */
   void parts(T& /* real */, T& /* imag */) const;

   /*
    * Complex conjugate.
    */
   friend complex<T,Syn> conj<T,Syn>(const complex<T,Syn>&);

protected:
   /*
    * Complex data.
    */
   T _real;       /* real component */
   T _imag;       /* imaginary component */
};

/***************************************************************************
 * Constructors and destructor.
 ***************************************************************************/

/*
 * Default constructor.
 * Return the complex number zero.
 */
template <typename T, typename Syn>
complex<T,Syn>::complex()
 : Syn(),
   _real(),
   _imag()
{ }

/*
 * Real constructor.
 * Convert a real number to complex form.
 */
template <typename T, typename Syn>
complex<T,Syn>::complex(const T& t)
 : Syn(),
   _real(t),
   _imag()
{ }

/*
 * Complex constructor.
 * Return a complex number with the given real and imaginary components.
 */
template <typename T, typename Syn>
complex<T,Syn>::complex(const T& t_real, const T& t_imag)
 : Syn(),
   _real(t_real),
   _imag(t_imag)
{ }

/*
 * Copy constructor.
 */
template <typename T, typename Syn>
complex<T,Syn>::complex(const complex<T,Syn>& z)
 : Syn() 
{
   auto_read_lock<const Syn> rlock(z);
   _real = z._real;
   _imag = z._imag;
}

/*
 * Type conversion copy constructor.
 */
template <typename T, typename Syn>
template <typename U, typename S>
complex<T,Syn>::complex(const complex<U,S>& z)
 : Syn() 
{
   auto_read_lock<const S> rlock(z);
   _real = z._real;
   _imag = z._imag;
}

/*
 * Destructor.
 */
template <typename T, typename Syn>
complex<T,Syn>::~complex() {
   /* do nothing */
}

/***************************************************************************
 * Named constructors.
 ***************************************************************************/

/*
 * Polar form.
 */
template <typename T, typename Syn>
complex<T,Syn> complex<T,Syn>::polar(const T& r, const T& theta) {
   const T r_copy(r);
   const T theta_copy(theta);
   return complex<T,Syn>(r_copy*cos(theta_copy), r_copy*sin(theta_copy));
}

/***************************************************************************
 * Serialization.
 ***************************************************************************/

/*
 * Serialize.
 */
template <typename T, typename Syn>
void complex<T,Syn>::serialize(
   serial_output_stream& s, const serializer<T>& slzr) const
{
   auto_read_lock<const Syn> rlock(*this);
   slzr.serialize(s, _real);
   slzr.serialize(s, _imag);
}

/*
 * Deserialize.
 */
template <typename T, typename Syn>
auto_ptr< complex<T,Syn> > complex<T,Syn>::deserialize(
   serial_input_stream& s, const serializer<T>& slzr) 
{
   auto_ptr<T> t_real = slzr.deserialize(s);
   auto_ptr<T> t_imag = slzr.deserialize(s);
   return auto_ptr< complex<T,Syn> >(new complex<T,Syn>(*t_real, *t_imag));
}

/***************************************************************************
 * I/O.
 ***************************************************************************/

/*
 * Formatted output to stream.
 */
template <typename T, typename Syn>
ostream& operator<<(ostream& os, const complex<T,Syn>& z) {
   auto_read_lock<const Syn> rlock(z);
   os << z._real;
   (z._imag < T()) ?
      (os << z._imag << 'i') : (os << '+' << abs(z._imag) << 'i'); 
   return os;
}

/***************************************************************************
 * Assignment operators: complex-real.
 ***************************************************************************/

template <typename T, typename Syn>
complex<T,Syn>& complex<T,Syn>::operator=(const T& t) {
   auto_write_lock<const Syn> wlock(*this);
   _real = t;
   _imag = T();
   return *this;
}

template <typename T, typename Syn>
complex<T,Syn>& complex<T,Syn>::operator+=(const T& t) {
   auto_write_lock<const Syn> wlock(*this);
   _real += t;
   return *this;
}

template <typename T, typename Syn>
complex<T,Syn>& complex<T,Syn>::operator-=(const T& t) {
   auto_write_lock<const Syn> wlock(*this);
   _real -= t;
   return *this;
}

template <typename T, typename Syn>
complex<T,Syn>& complex<T,Syn>::operator*=(const T& t) {
   const T t_copy(t);
   auto_write_lock<const Syn> wlock(*this);
   _real *= t_copy;
   _imag *= t_copy;
   return *this;
}

template <typename T, typename Syn>
complex<T,Syn>& complex<T,Syn>::operator/=(const T& t) {
   const T t_copy(t);
   auto_write_lock<const Syn> wlock(*this);
   _real /= t_copy;
   _imag /= t_copy;
   return *this;
}

/***************************************************************************
 * Assignment operators: complex-complex.
 ***************************************************************************/

template <typename T, typename Syn>
complex<T,Syn>& complex<T,Syn>::operator=(const complex<T,Syn>& z) {
   auto_read_write_lock<const Syn> rwlock(z, *this);
   _real = z._real;
   _imag = z._imag;
   return *this;
}

template <typename T, typename Syn>
complex<T,Syn>& complex<T,Syn>::operator+=(const complex<T,Syn>& z) {
   auto_read_write_lock<const Syn> rwlock(z, *this);
   _real += z._real;
   _imag += z._imag;
   return *this;
}

template <typename T, typename Syn>
complex<T,Syn>& complex<T,Syn>::operator-=(const complex<T,Syn>& z) {
   auto_read_write_lock<const Syn> rwlock(z, *this);
   _real -= z._real;
   _imag -= z._imag;
   return *this;
}

template <typename T, typename Syn>
complex<T,Syn>& complex<T,Syn>::operator*=(const complex<T,Syn>& z) {
   auto_read_write_lock<const Syn> rwlock(z, *this);
   const T this_real(_real);
   const T this_imag(_imag);
   _real = (this_real*z._real) - (this_imag*z._imag);
   _imag = (this_real*z._imag) + (this_imag*z._real);
   return *this;
}

template <typename T, typename Syn>
complex<T,Syn>& complex<T,Syn>::operator/=(const complex<T,Syn>& z) {
   auto_read_write_lock<const Syn> rwlock(z, *this);
   const T this_real(_real);
   const T this_imag(_imag);
   T mag = (z._real*z._real) + (z._imag*z._imag);
   _real = ((this_real*z._real) + (this_imag*z._imag))/mag;
   _imag = ((this_imag*z._real) - (this_real*z._imag))/mag;
   return *this;
}

template <typename T, typename Syn>
template <typename U, typename S>
complex<T,Syn>& complex<T,Syn>::operator=(const complex<U,S>& z) {
   auto_read_write_lock<const synchronizable> rwlock(z, *this);
   _real = z._real;
   _imag = z._imag;
   return *this;
}

template <typename T, typename Syn>
template <typename U, typename S>
complex<T,Syn>& complex<T,Syn>::operator+=(const complex<U,S>& z) {
   auto_read_write_lock<const synchronizable> rwlock(z, *this);
   _real += z._real;
   _imag += z._imag;
   return *this;
}

template <typename T, typename Syn>
template <typename U, typename S>
complex<T,Syn>& complex<T,Syn>::operator-=(const complex<U,S>& z) {
   auto_read_write_lock<const synchronizable> rwlock(z, *this);
   _real -= z._real;
   _imag -= z._imag;
   return *this;
}

template <typename T, typename Syn>
template <typename U, typename S>
complex<T,Syn>& complex<T,Syn>::operator*=(const complex<U,S>& z) {
   auto_read_write_lock<const synchronizable> rwlock(z, *this);
   const T this_real(_real);
   const T this_imag(_imag);
   _real = (this_real*z._real) - (this_imag*z._imag);
   _imag = (this_real*z._imag) + (this_imag*z._real);
   return *this;
}

template <typename T, typename Syn>
template <typename U, typename S>
complex<T,Syn>& complex<T,Syn>::operator/=(const complex<U,S>& z) {
   auto_read_write_lock<const synchronizable> rwlock(z, *this);
   const T this_real(_real);
   const T this_imag(_imag);
   U mag = (z._real*z._real) + (z._imag*z._imag);
   /* must divide by mag below (not beforehand) */
   /* because of type promotion rules           */
   _real = ((this_real*z._real) + (this_imag*z._imag))/mag;
   _imag = ((this_imag*z._real) - (this_real*z._imag))/mag;
   return *this;
}

/***************************************************************************
 * Binary complex-real operators.
 ***************************************************************************/

template <typename T, typename Syn>
inline complex<T,Syn> operator+(const complex<T,Syn>& z, const T& t) {
   return (complex<T,Syn>(z) += t);
}

template <typename T, typename Syn>
inline complex<T,Syn> operator-(const complex<T,Syn>& z, const T& t) {
   return (complex<T,Syn>(z) -= t);
}

template <typename T, typename Syn>
inline complex<T,Syn> operator*(const complex<T,Syn>& z, const T& t) {
   return (complex<T,Syn>(z) *= t);
}  

template <typename T, typename Syn>
inline complex<T,Syn> operator/(const complex<T,Syn>& z, const T& t) {
   return (complex<T,Syn>(z) /= t);
}  

template <typename T, typename Syn>
inline complex<T,Syn> operator+(const T& t, const complex<T,Syn>& z) {
   return (complex<T,Syn>(z) += t);
}
   
template <typename T, typename Syn>
inline complex<T,Syn> operator-(const T& t, const complex<T,Syn>& z) {
   return (complex<T,Syn>(t) -= z);
}

template <typename T, typename Syn>
inline complex<T,Syn> operator*(const T& t, const complex<T,Syn>& z) {
   return (complex<T,Syn>(z) *= t);
}  

template <typename T, typename Syn>
inline complex<T,Syn> operator/(const T& t, const complex<T,Syn>& z) {
   return (complex<T,Syn>(t) /= z); 
}  

/***************************************************************************
 * Binary complex-complex operators.
 ***************************************************************************/

template <typename T, typename Syn>
inline complex<T,Syn> operator+(
   const complex<T,Syn>& z0, const complex<T,Syn>& z1)
{
   return (complex<T,Syn>(z0) += z1);
}

template <typename T, typename Syn>
inline complex<T,Syn> operator-(
   const complex<T,Syn>& z0, const complex<T,Syn>& z1)
{
   return (complex<T,Syn>(z0) -= z1);
}

template <typename T, typename Syn>
inline complex<T,Syn> operator*(
   const complex<T,Syn>& z0, const complex<T,Syn>& z1)
{
   return (complex<T,Syn>(z0) *= z1);
}

template <typename T, typename Syn>
inline complex<T,Syn> operator/(
   const complex<T,Syn>& z0, const complex<T,Syn>& z1)
{
   return (complex<T,Syn>(z0) /= z1);
}
 
/***************************************************************************
 * Binary complex-real comparators.
 *
 * The == and != operators return whether the complex number is/is not 
 * exactly the given real number (with zero imaginary component).  
 *
 * The other operators (less than/greater than judgements) compare the real
 * part of the complex number with the given real number and only use the 
 * imaginary component in the event of a tie.  A complex number z is greater 
 * than a real number t iff (real(z) > t) or (real(z) == t and imag(z) > 0).
 ***************************************************************************/

template <typename T, typename Syn>
bool operator==(const complex<T,Syn>& z, const T& t) {
   auto_read_lock<const Syn> rlock(z);
   return ((z._real == t) && (z._imag == T()));
}

template <typename T, typename Syn>
bool operator!=(const complex<T,Syn>& z, const T& t) {
   auto_read_lock<const Syn> rlock(z);
   return ((z._real != t) || (z._imag != T()));
}

template <typename T, typename Syn>
bool operator<(const complex<T,Syn>& z, const T& t) {
   const T t_copy(t);
   auto_read_lock<const Syn> rlock(z);
   return ((z._real < t_copy) || ((z._real == t_copy) && (z._imag < T())));
}

template <typename T, typename Syn>
bool operator>(const complex<T,Syn>& z, const T& t) {
   const T t_copy(t);
   auto_read_lock<const Syn> rlock(z);
   return ((z._real > t_copy) || ((z._real == t_copy) && (z._imag > T())));
}

template <typename T, typename Syn>
bool operator<=(const complex<T,Syn>& z, const T& t) {
   const T t_copy(t);
   auto_read_lock<const Syn> rlock(z);
   return ((z._real < t_copy) || ((z._real == t_copy) && (z._imag <= T())));
}

template <typename T, typename Syn>
bool operator>=(const complex<T,Syn>& z, const T& t) {
   const T t_copy(t);
   auto_read_lock<const Syn> rlock(z);
   return ((z._real > t_copy) || ((z._real == t_copy) && (z._imag >= T())));
}

template <typename T, typename Syn>
inline bool operator==(const T& t, const complex<T,Syn>& z) {
   return (z == t);
}

template <typename T, typename Syn>
inline bool operator!=(const T& t, const complex<T,Syn>& z) {
   return (z != t);
}

template <typename T, typename Syn>
inline bool operator<(const T& t, const complex<T,Syn>& z) {
   return (z > t);
}

template <typename T, typename Syn>
inline bool operator>(const T& t, const complex<T,Syn>& z) {
   return (z < t);
}

template <typename T, typename Syn>
inline bool operator<=(const T& t, const complex<T,Syn>& z) {
   return (z >= t);
}

template <typename T, typename Syn>
inline bool operator>=(const T& t, const complex<T,Syn>& z) {
   return (z <= t);
}

/***************************************************************************
 * Binary complex-complex comparators.
 *
 * The == and != operators return whether the two complex numbers are/are not 
 * exactly equal (both real and imaginary components equal/not equal).
 *
 * The other operators (less than/greater than judgements) compare the real 
 * parts of the two complex numbers and use the imaginary components in the 
 * event of a tie.  Complex number z0 is greater than complex number z1 iff 
 * (real(z0) > real(z1)) or (real(z0) == real(z1) and imag(z0) > imag(z1)).
 ***************************************************************************/

template <typename T, typename Syn>
bool operator==(const complex<T,Syn>& z0, const complex<T,Syn>& z1) {
   auto_read_read_lock<const Syn> rrlock(z0, z1);
   return (((z0._real == z1._real) && (z0._imag == z1._imag)));
}

template <typename T, typename Syn>
bool operator!=(const complex<T,Syn>& z0, const complex<T,Syn>& z1) {
   auto_read_read_lock<const Syn> rrlock(z0, z1);
   return (((z0._real != z1._real) || (z0._imag != z1._imag)));
}

template <typename T, typename Syn>
bool operator<(const complex<T,Syn>& z0, const complex<T,Syn>& z1) {
   auto_read_read_lock<const Syn> rrlock(z0, z1);
   return (
      (z0._real < z1._real) ||
      ((z0._real == z1._real) && (z0._imag < z1._imag))
   );
}

template <typename T, typename Syn>
bool operator>(const complex<T,Syn>& z0, const complex<T,Syn>& z1) {
   auto_read_read_lock<const Syn> rrlock(z0, z1);
   return (
      (z0._real > z1._real) ||
      ((z0._real == z1._real) && (z0._imag > z1._imag))
   );
}

template <typename T, typename Syn>
bool operator<=(const complex<T,Syn>& z0, const complex<T,Syn>& z1) {
   auto_read_read_lock<const Syn> rrlock(z0, z1);
   return (
      (z0._real < z1._real) ||
      ((z0._real == z1._real) && (z0._imag <= z1._imag))
   );
}

template <typename T, typename Syn>
bool operator>=(const complex<T,Syn>& z0, const complex<T,Syn>& z1) {
   auto_read_read_lock<const Syn> rrlock(z0, z1);
   return (
      (z0._real > z1._real) ||
      ((z0._real == z1._real) && (z0._imag >= z1._imag))
   );
}

/***************************************************************************
 * Unary operators.
 ***************************************************************************/

/*
 * Unary arithmetic operators.
 */
template <typename T, typename Syn>
complex<T,Syn> complex<T,Syn>::operator+() const {
   auto_read_lock<const Syn> rlock(*this);
   return complex<T,Syn>(+(_real), +(_imag));
}
   
template <typename T, typename Syn>
complex<T,Syn> complex<T,Syn>::operator-() const {
   auto_read_lock<const Syn> rlock(*this);
   return complex<T,Syn>(-(_real), -(_imag));
}

template <typename T, typename Syn>
complex<T,Syn>& complex<T,Syn>::operator++() {
   auto_write_lock<const Syn> wlock(*this);
   ++(_real);
   return *this;
}

template <typename T, typename Syn>
complex<T,Syn>& complex<T,Syn>::operator--() {
   auto_write_lock<const Syn> wlock(*this);
   --(_real);
   return *this;
}

template <typename T, typename Syn>
complex<T,Syn> complex<T,Syn>::operator++(int) {
   auto_write_lock<const Syn> wlock(*this);
   return complex<T,Syn>(_real++, _imag);
}

template <typename T, typename Syn>
complex<T,Syn> complex<T,Syn>::operator--(int) {
   auto_write_lock<const Syn> wlock(*this);
   return complex<T,Syn>(_real--, _imag);
}

/*
 * Unary logical inversion operator.
 */
template <typename T, typename Syn>
inline bool complex<T,Syn>::operator!() const {
   return (*this == T());
}

/*
 * Unary multiplicative inversion operator.
 */
template <typename T, typename Syn>
complex<T,Syn> complex<T,Syn>::operator~() const {
   auto_read_lock<const Syn> rlock(*this);
   T mag = (_real*_real) + (_imag*_imag);
   return complex<T,Syn>((_real/mag), (-_imag/mag));
}

/***************************************************************************
 * Transcendentals.
 ***************************************************************************/

/*
 * Cosine.
 */
template <typename T, typename Syn>
complex<T,Syn> cos(const complex<T,Syn>& z) {
   auto_read_lock<const Syn> rlock(z);
   return complex<T,Syn>(
      cos(z._real)*cosh(z._imag), -sin(z._real)*sinh(z._imag)
   );
}

/*
 * Sine.
 */
template <typename T, typename Syn>
complex<T,Syn> sin(const complex<T,Syn>& z) {
   auto_read_lock<const Syn> rlock(z);
   return complex<T,Syn>(
      sin(z._real)*cosh(z._imag), cos(z._real)*sinh(z._imag)
   );
}

/*
 * Tangent.
 */
template <typename T, typename Syn>
complex<T,Syn> tan(const complex<T,Syn>& z) {
   auto_read_lock<const Syn> rlock(z);
   T cos_x = cos(z._real);
   T sin_x = sin(z._real);
   T cosh_y = cosh(z._imag);
   T sinh_y = sinh(z._imag);
   T x = sin_x*cosh_y;
   T y = cos_x*sinh_y;
   T a = cos_x*cosh_y;
   T b = -sin_x*sinh_y;
   T mag = a*a + b*b;
   return complex<T,Syn>((x*a + y*b)/mag, (y*a - x*b)/mag);
}

/*
 * Hyperbolic cosine.
 */
template <typename T, typename Syn>
complex<T,Syn> cosh(const complex<T,Syn>& z) {
   auto_read_lock<const Syn> rlock(z);
   return complex<T,Syn>(
      cosh(z._real)*cos(z._imag), sinh(z._real)*sin(z._imag)
   );
}

/*
 * Hyperbolic sine.
 */
template <typename T, typename Syn>
complex<T,Syn> sinh(const complex<T,Syn>& z) {
   auto_read_lock<const Syn> rlock(z);
   return complex<T,Syn>(
      sinh(z._real)*cos(z._imag), cosh(z._real)*sin(z._imag)
   );
}

/*
 * Hyperbolic tangent.
 */
template <typename T, typename Syn>
complex<T,Syn> tanh(const complex<T,Syn>& z) {
   auto_read_lock<const Syn> rlock(z);
   T cosh_x = cosh(z._real);
   T sinh_x = sinh(z._real);
   T cos_y = cos(z._imag);
   T sin_y = sin(z._imag);
   T x = sinh_x*cos_y;
   T y = cosh_x*sin_y;
   T a = cosh_x*cos_y;
   T b = sinh_x*sin_y;
   T mag = a*a + b*b;
   return complex<T,Syn>((x*a + y*b)/mag, (y*a - x*b)/mag);
}

/*
 * Arccosine (principal value).
 */
template <typename T, typename Syn>
complex<T,Syn> acos(const complex<T,Syn>& z) {
   const complex<T,Syn> z_copy(z);
   return (-(T(M_PIl)) * log(z_copy + T(M_PIl)*sqrt(T(1) - z_copy*z_copy)));
}

/*
 * Arcsine (principal value).
 */
template <typename T, typename Syn>
complex<T,Syn> asin(const complex<T,Syn>& z) {
   const complex<T,Syn> z_copy(z);
   return (-(T(M_PIl)) * log(T(M_PIl)*z_copy + sqrt(T(1) - z_copy*z_copy)));
}

/*
 * Arctangent (principal value).
 */
template <typename T, typename Syn>
complex<T,Syn> atan(const complex<T,Syn>& z) {
   const complex<T,Syn> z_copy(z);
   return (T(M_PI_2l) * log((T(M_PIl) + z_copy) / (T(M_PIl) - z_copy)));
}

/*
 * Hyperbolic arccosine (principal value).
 */
template <typename T, typename Syn>
complex<T,Syn> acosh(const complex<T,Syn>& z) {
   const complex<T,Syn> z_copy(z);
   return log(z_copy + sqrt(z_copy*z_copy - T(1)));
}

/*
 * Hyperbolic arcsine (principal value).
 */
template <typename T, typename Syn>
complex<T,Syn> asinh(const complex<T,Syn>& z) {
   const complex<T,Syn> z_copy(z);
   return log(z_copy + sqrt(z_copy*z_copy + T(1)));
}

/*
 * Hyperbolic arctangent (principal value).
 */
template <typename T, typename Syn>
complex<T,Syn> atanh(const complex<T,Syn>& z) {
   const complex<T,Syn> z_copy(z);
   return (log((T(1) + z_copy)/(T(1) - z_copy))/T(2));
}

/*
 * Square root (principal value).
 */
template <typename T, typename Syn>
complex<T,Syn> sqrt(const complex<T,Syn>& z) {
   /* read source */
   auto_read_lock<const Syn> rlock(z);
   const T x(z._real);
   const T y(z._imag);
   /* compute square root */
   if (x == T()) {
      T t = sqrt(abs(y)/T(2));
      return complex<T,Syn>(t, (y < T()) ? (-t) : t);
   } else {
      T t = sqrt(T(2) * (sqrt(x*x + y*y) + abs(x)));
      T u = t / T(2);
      return (x > T()) ? complex<T,Syn>(u, y/t) :
                         complex<T,Syn>(abs(y)/t, (y < T()) ? (-u) : u);
   }
}

/*
 * Log (base e) (principal value).
 */
template <typename T, typename Syn>
complex<T,Syn> log(const complex<T,Syn>& z) {
   const complex<T,Syn> z_copy(z);
   return complex<T,Syn>(log(abs(z_copy)), arg(z_copy));
}

/*
 * Exponential.
 */
template <typename T, typename Syn>
complex<T,Syn> exp(const complex<T,Syn>& z) {
   auto_read_lock<const Syn> rlock(z);
   T e_x = exp(z._real);
   return complex<T,Syn>(e_x*cos(z._imag), e_x*sin(z._imag));
}

/*
 * Power (real raised to complex power).
 */
template <typename T, typename Syn>
complex<T,Syn> pow(const T& t, const complex<T,Syn>& z) {
   return exp(z * log(t));
}

/*
 * Power (integer).
 */
template <typename T, typename Syn>
complex<T,Syn> pow(const complex<T,Syn>& z, int n) {
   return exp(T(n) * log(z));
}

/*
 * Power (real).
 */
template <typename T, typename Syn>
complex<T,Syn> pow(const complex<T,Syn>& z, const T& t) {
   return exp(t * log(z));
}

/*
 * Power (complex).
 */
template <typename T, typename Syn>
complex<T,Syn> pow(const complex<T,Syn>& z_base, const complex<T,Syn>& z_exp) {
   return exp(z_exp * log(z_base));
}

/*
 * Two-parameter arctangent of the real parts.
 */
template <typename T, typename Syn>
complex<T,Syn> atan2(const complex<T,Syn>& zy, const complex<T,Syn>& zx) {
   return atan2(zy.real(), zx.real());
}

/***************************************************************************
 * Polar form.
 ***************************************************************************/

/*
 * Magnitude.
 */
template <typename T, typename Syn>
T abs(const complex<T,Syn>& z) {
   auto_read_lock<const Syn> rlock(z);
   return sqrt((z._real*z._real) + (z._imag*z._imag));
}

/*
 * Argument.
 */
template <typename T, typename Syn>
T arg(const complex<T,Syn>& z) {
   auto_read_lock<const Syn> rlock(z);
   return atan2(z._imag, z._real);
}

/***************************************************************************
 * Real/imaginary components.
 ***************************************************************************/

/*
 * Check if real-valued (zero imaginary component).
 */
template <typename T, typename Syn>
bool complex<T,Syn>::is_real() const {
   auto_read_lock<const Syn> rlock(*this);
   return (_imag == T());
}

/*
 * Check if imaginary-valued (zero real component). 
 */
template <typename T, typename Syn>
bool complex<T,Syn>::is_imag() const {
   auto_read_lock<const Syn> rlock(*this);
   return (_real == T());
}

/*
 * Real component.
 */
template <typename T, typename Syn>
T complex<T,Syn>::real() const {
   auto_read_lock<const Syn> rlock(*this);
   return _real;
}

/*
 * Imaginary component. 
 */
template <typename T, typename Syn>
T complex<T,Syn>::imag() const {
   auto_read_lock<const Syn> rlock(*this);
   return _imag;
}

/*
 * Real and imaginary components.
 */
template <typename T, typename Syn>
void complex<T,Syn>::parts(T& t_real, T& t_imag) const {
   auto_read_lock<const Syn> rlock(*this);
   t_real = _real;
   t_imag = _imag;
}

/*
 * Complex conjugate.
 */
template <typename T, typename Syn>
complex<T,Syn> conj(const complex<T,Syn>& z) {
   auto_read_lock<const Syn> rlock(z);
   return complex<T,Syn>(z._real, -(z._imag));
}

} /* namespace math */

#endif
