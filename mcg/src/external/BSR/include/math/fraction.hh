/*
 * Fractions (thread-safe).
 */
#ifndef MATH__FRACTION_HH
#define MATH__FRACTION_HH

#include "concurrent/threads/synchronization/locks/auto_read_lock.hh"
#include "concurrent/threads/synchronization/locks/auto_read_read_lock.hh"
#include "concurrent/threads/synchronization/locks/auto_read_write_lock.hh"
#include "concurrent/threads/synchronization/locks/auto_write_lock.hh"
#include "concurrent/threads/synchronization/synchronizables/synchronizable.hh"
#include "concurrent/threads/synchronization/synchronizables/unsynchronized.hh"
#include "interfaces/comparable.hh"
#include "io/serialization/serial_input_stream.hh"
#include "io/serialization/serial_output_stream.hh"
#include "io/serialization/serializers.hh"
#include "io/streams/ostream.hh"
#include "lang/pointers/auto_ptr.hh"

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
using interfaces::comparable;
using io::serialization::serial_input_stream;
using io::serialization::serial_output_stream;
using io::serialization::serializer;
using io::serialization::serializers;
using io::streams::ostream;
using lang::pointers::auto_ptr;

/*
 * Declare existence of fraction class.
 */
template <typename T, typename Syn> class fraction;

/*
 * Declare prototypes for fraction template friend functions.
 */
template <typename T, typename Syn>
ostream& operator<<(ostream&, const fraction<T,Syn>&);

template <typename T, typename Syn>
fraction<T,Syn> operator+(const fraction<T,Syn>&, const T&);
template <typename T, typename Syn>
fraction<T,Syn> operator-(const fraction<T,Syn>&, const T&);
template <typename T, typename Syn>
fraction<T,Syn> operator*(const fraction<T,Syn>&, const T&);
template <typename T, typename Syn>
fraction<T,Syn> operator/(const fraction<T,Syn>&, const T&);

template <typename T, typename Syn>
fraction<T,Syn> operator+(const T&, const fraction<T,Syn>&);
template <typename T, typename Syn>
fraction<T,Syn> operator-(const T&, const fraction<T,Syn>&);
template <typename T, typename Syn>
fraction<T,Syn> operator*(const T&, const fraction<T,Syn>&);
template <typename T, typename Syn>
fraction<T,Syn> operator/(const T&, const fraction<T,Syn>&);

template <typename T, typename Syn>
fraction<T,Syn> operator+(const fraction<T,Syn>&, const fraction<T,Syn>&);
template <typename T, typename Syn>
fraction<T,Syn> operator-(const fraction<T,Syn>&, const fraction<T,Syn>&);
template <typename T, typename Syn>
fraction<T,Syn> operator*(const fraction<T,Syn>&, const fraction<T,Syn>&);
template <typename T, typename Syn>
fraction<T,Syn> operator/(const fraction<T,Syn>&, const fraction<T,Syn>&);

template <typename T, typename Syn>
bool operator==(const fraction<T,Syn>&, const T&);
template <typename T, typename Syn>
bool operator!=(const fraction<T,Syn>&, const T&);
template <typename T, typename Syn>
bool operator< (const fraction<T,Syn>&, const T&);
template <typename T, typename Syn>
bool operator> (const fraction<T,Syn>&, const T&);
template <typename T, typename Syn>
bool operator<=(const fraction<T,Syn>&, const T&);
template <typename T, typename Syn>
bool operator>=(const fraction<T,Syn>&, const T&);

template <typename T, typename Syn>
bool operator==(const T&, const fraction<T,Syn>&);
template <typename T, typename Syn>
bool operator!=(const T&, const fraction<T,Syn>&);
template <typename T, typename Syn>
bool operator< (const T&, const fraction<T,Syn>&);
template <typename T, typename Syn>
bool operator> (const T&, const fraction<T,Syn>&);
template <typename T, typename Syn>
bool operator<=(const T&, const fraction<T,Syn>&);
template <typename T, typename Syn>
bool operator>=(const T&, const fraction<T,Syn>&);

template <typename T, typename Syn>
bool operator==(const fraction<T,Syn>&, const fraction<T,Syn>&);
template <typename T, typename Syn>
bool operator!=(const fraction<T,Syn>&, const fraction<T,Syn>&);
template <typename T, typename Syn>
bool operator< (const fraction<T,Syn>&, const fraction<T,Syn>&);
template <typename T, typename Syn>
bool operator> (const fraction<T,Syn>&, const fraction<T,Syn>&);
template <typename T, typename Syn>
bool operator<=(const fraction<T,Syn>&, const fraction<T,Syn>&);
template <typename T, typename Syn>
bool operator>=(const fraction<T,Syn>&, const fraction<T,Syn>&);

template <typename T, typename Syn>
fraction<T,Syn> abs(const fraction<T,Syn>&);

/*
 * Fraction class.
 */
template <typename T, typename Syn = unsynchronized>
class fraction : public comparable< fraction<T,Syn> >,
                 protected Syn {
public:
   /*
    * Friend classes.
    */
   template <typename U, typename S> friend class fraction;

   /*
    * Constructors.
    */
   fraction();
   explicit fraction(const T&);
   explicit fraction(const T&, const T&);

   /*
    * Copy constructors.
    */
   fraction(const fraction<T,Syn>&);
   template <typename U, typename S> fraction(const fraction<U,S>&);
    
   /*
    * Destructor.
    */
   virtual ~fraction();

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
   static auto_ptr< fraction<T,Syn> > deserialize(
      serial_input_stream&,
      const serializer<T>& = serializers<T>::s_default()
   );

   /*
    * Formatted output to stream.
    */
   friend ostream& operator<< <T,Syn>(ostream&, const fraction<T,Syn>&);
   
   /*
    * Assignment operators: fraction-value.
    */
   fraction<T,Syn>& operator=(const T&);
   fraction<T,Syn>& operator+=(const T&);
   fraction<T,Syn>& operator-=(const T&);
   fraction<T,Syn>& operator*=(const T&);
   fraction<T,Syn>& operator/=(const T&);
  
   /*
    * Assignment operators: fraction-fraction.
    */
   fraction<T,Syn>& operator=(const fraction<T,Syn>&);
   fraction<T,Syn>& operator+=(const fraction<T,Syn>&);
   fraction<T,Syn>& operator-=(const fraction<T,Syn>&);
   fraction<T,Syn>& operator*=(const fraction<T,Syn>&);
   fraction<T,Syn>& operator/=(const fraction<T,Syn>&);
   
   template <typename U, typename S> fraction<T,Syn>&
      operator=(const fraction<U,S>&);
   template <typename U, typename S> fraction<T,Syn>&
      operator+=(const fraction<U,S>&);
   template <typename U, typename S> fraction<T,Syn>&
      operator-=(const fraction<U,S>&);
   template <typename U, typename S> fraction<T,Syn>&
      operator*=(const fraction<U,S>&);
   template <typename U, typename S> fraction<T,Syn>&
      operator/=(const fraction<U,S>&);
   
   /*
    * Binary fraction-value operators.
    */
   friend fraction<T,Syn> operator+ <T,Syn>(const fraction<T,Syn>&, const T&);
   friend fraction<T,Syn> operator- <T,Syn>(const fraction<T,Syn>&, const T&);
   friend fraction<T,Syn> operator* <T,Syn>(const fraction<T,Syn>&, const T&);
   friend fraction<T,Syn> operator/ <T,Syn>(const fraction<T,Syn>&, const T&);

   friend fraction<T,Syn> operator+ <T,Syn>(const T&, const fraction<T,Syn>&);
   friend fraction<T,Syn> operator- <T,Syn>(const T&, const fraction<T,Syn>&);
   friend fraction<T,Syn> operator* <T,Syn>(const T&, const fraction<T,Syn>&);
   friend fraction<T,Syn> operator/ <T,Syn>(const T&, const fraction<T,Syn>&);
   
   /*
    * Binary fraction-fraction operators.
    */
   friend fraction<T,Syn>
      operator+ <T,Syn>(const fraction<T,Syn>&, const fraction<T,Syn>&);
   friend fraction<T,Syn>
      operator- <T,Syn>(const fraction<T,Syn>&, const fraction<T,Syn>&);
   friend fraction<T,Syn>
      operator* <T,Syn>(const fraction<T,Syn>&, const fraction<T,Syn>&);
   friend fraction<T,Syn>
      operator/ <T,Syn>(const fraction<T,Syn>&, const fraction<T,Syn>&);

   /*
    * Binary fraction-value comparators.
    */
   friend bool operator== <T,Syn>(const fraction<T,Syn>&, const T&);
   friend bool operator!= <T,Syn>(const fraction<T,Syn>&, const T&);
   friend bool operator<  <T,Syn>(const fraction<T,Syn>&, const T&);
   friend bool operator>  <T,Syn>(const fraction<T,Syn>&, const T&);
   friend bool operator<= <T,Syn>(const fraction<T,Syn>&, const T&);
   friend bool operator>= <T,Syn>(const fraction<T,Syn>&, const T&);
   
   friend bool operator== <T,Syn>(const T&, const fraction<T,Syn>&);
   friend bool operator!= <T,Syn>(const T&, const fraction<T,Syn>&);
   friend bool operator<  <T,Syn>(const T&, const fraction<T,Syn>&);
   friend bool operator>  <T,Syn>(const T&, const fraction<T,Syn>&);
   friend bool operator<= <T,Syn>(const T&, const fraction<T,Syn>&);
   friend bool operator>= <T,Syn>(const T&, const fraction<T,Syn>&);
   
   /*
    * Binary fraction-fraction comparators.
    */
   friend bool
      operator== <T,Syn>(const fraction<T,Syn>&, const fraction<T,Syn>&);
   friend bool
      operator!= <T,Syn>(const fraction<T,Syn>&, const fraction<T,Syn>&);
   friend bool
      operator<  <T,Syn>(const fraction<T,Syn>&, const fraction<T,Syn>&);
   friend bool
      operator>  <T,Syn>(const fraction<T,Syn>&, const fraction<T,Syn>&);
   friend bool
      operator<= <T,Syn>(const fraction<T,Syn>&, const fraction<T,Syn>&);
   friend bool
      operator>= <T,Syn>(const fraction<T,Syn>&, const fraction<T,Syn>&);

   /*
    * Comparison method.
    */
   int compare_to(const fraction<T,Syn>&) const;

   /*
    * Unary arithmetic operators.
    */
   fraction<T,Syn> operator+() const;
   fraction<T,Syn> operator-() const;
   fraction<T,Syn>& operator++();
   fraction<T,Syn>& operator--();
   fraction<T,Syn> operator++(int);
   fraction<T,Syn> operator--(int);

   /*
    * Unary logical inversion operator.
    */
   bool operator!() const;

   /*
    * Unary multiplicative inversion operator.
    */
   fraction<T,Syn> operator~() const;

   /*
    * Return the sign (-1, 0, or 1) of the value represented by the fraction.
    */
   int sign() const;

   /*
    * Absolute value.
    */
   friend fraction<T,Syn> abs<T,Syn>(const fraction<T,Syn>&);

   /*
    * Numerator.
    */
   T numerator() const;

   /*
    * Denominator.
    */
   T denominator() const;

   /*
    * Value.
    * Return the numerator divided by the denominator.
    */
   T value() const;

   /*
    * Reduce the fraction.
    * If T is an integer type, divide both the numerator and denominator by
    * their greatest common divisor.  Otherwise, if the denominator is nonzero,
    * divide the numerator by the denominator and set the denominator to one.
    */
   fraction<T,Syn>& reduce();

protected:
   /*
    * Fraction data.
    */
   T _numerator;     /* numerator */
   T _denominator;   /* denominator */
};

/***************************************************************************
 * Constructors and destructor.
 ***************************************************************************/

/*
 * Constructor.
 * Return the fraction zero.
 */
template <typename T, typename Syn>
fraction<T,Syn>::fraction()
 : Syn(),
   _numerator(),
   _denominator(1)
{ }

/*
 * Constructor.
 * Return the given value.
 */
template <typename T, typename Syn>
fraction<T,Syn>::fraction(const T& t)
 : Syn(),
   _numerator(t),
   _denominator(1)
{ }

/*
 * Constructor.
 * Return a fraction with the given numerator and denominator.
 */
template <typename T, typename Syn>
fraction<T,Syn>::fraction(const T& num, const T& denom)
 : Syn(),
   _numerator(num),
   _denominator(denom)
{ }

/*
 * Copy constructor.
 */
template <typename T, typename Syn>
fraction<T,Syn>::fraction(const fraction<T,Syn>& f)
 : Syn(),
   _numerator(f._numerator),
   _denominator(f._denominator)
{ }

/*
 * Type conversion copy constructor.
 */
template <typename T, typename Syn>
template <typename U, typename S>
fraction<T,Syn>::fraction(const fraction<U,S>& f)
 : Syn()
{
   auto_read_lock<const S> rlock(f);
   _numerator   = f._numerator;
   _denominator = f._denominator;
}
 
/*
 * Destructor.
 */
template <typename T, typename Syn>
fraction<T,Syn>::~fraction() {
   /* do nothing */
}

/***************************************************************************
 * Serialization.
 ***************************************************************************/

/*
 * Serialize.
 */
template <typename T, typename Syn>
void fraction<T,Syn>::serialize(
   serial_output_stream& s, const serializer<T>& slzr) const
{
   auto_read_lock<const Syn> rlock(*this);
   slzr.serialize(s, _numerator);
   slzr.serialize(s, _denominator);
}

/*
 * Deserialize.
 */
template <typename T, typename Syn>
auto_ptr< fraction<T,Syn> > fraction<T,Syn>::deserialize(
  serial_input_stream& s, const serializer<T>& slzr) 
{
   auto_ptr<T> num   = slzr.deserialize(s);
   auto_ptr<T> denom = slzr.deserialize(s);
   return auto_ptr< fraction<T,Syn> >(new fraction<T,Syn>(*num, *denom));
}

/***************************************************************************
 * I/O.
 ***************************************************************************/

/*
 * Formatted output to stream.
 */
template <typename T, typename Syn>
ostream& operator<<(ostream& os, const fraction<T,Syn>& f) {
   auto_read_lock<const Syn> rlock(f);
   os << f._numerator;
   if (f._denominator != T(1))
      os << '/' << f._denominator;
   return os;
}

/***************************************************************************
 * Assignment operators: fraction-value.
 ***************************************************************************/

template <typename T, typename Syn>
fraction<T,Syn>& fraction<T,Syn>::operator=(const T& t) {
   auto_write_lock<const Syn> wlock(*this);
   _numerator   = t;
   _denominator = T(1);
   return *this;
}

template <typename T, typename Syn>
fraction<T,Syn>& fraction<T,Syn>::operator+=(const T& t) {
   auto_write_lock<const Syn> wlock(*this);
   _numerator += (t * _denominator);
   return *this;
}

template <typename T, typename Syn>
fraction<T,Syn>& fraction<T,Syn>::operator-=(const T& t) {
   auto_write_lock<const Syn> wlock(*this);
   _numerator -= (t * _denominator);
   return *this;
}

template <typename T, typename Syn>
fraction<T,Syn>& fraction<T,Syn>::operator*=(const T& t) {
   auto_write_lock<const Syn> wlock(*this);
   _numerator *= t;
   return *this;
}

template <typename T, typename Syn>
fraction<T,Syn>& fraction<T,Syn>::operator/=(const T& t) {
   auto_write_lock<const Syn> wlock(*this);
   _denominator *= t;
   return *this;
}

/***************************************************************************
 * Assignment operators: fraction-fraction.
 ***************************************************************************/

template <typename T, typename Syn>
fraction<T,Syn>& fraction<T,Syn>::operator=(const fraction<T,Syn>& f) {
   auto_read_write_lock<const Syn> rwlock(f, *this);
   _numerator   = f._numerator;
   _denominator = f._denominator;
   return *this;
}

template <typename T, typename Syn>
fraction<T,Syn>& fraction<T,Syn>::operator+=(const fraction<T,Syn>& f) {
   auto_read_write_lock<const Syn> rwlock(f, *this);
   _numerator = (_numerator * f._denominator) + (f._numerator * _denominator);
   _denominator *= f._denominator;
   return *this;
}

template <typename T, typename Syn>
fraction<T,Syn>& fraction<T,Syn>::operator-=(const fraction<T,Syn>& f) {
   auto_read_write_lock<const Syn> rwlock(f, *this);
   _numerator = (_numerator * f._denominator) - (f._numerator * _denominator);
   _denominator *= f._denominator;
   return *this;
}

template <typename T, typename Syn>
fraction<T,Syn>& fraction<T,Syn>::operator*=(const fraction<T,Syn>& f) {
   auto_read_write_lock<const Syn> rwlock(f, *this);
   _numerator   *= f._numerator;
   _denominator *= f._denominator;
   return *this;
}

template <typename T, typename Syn>
fraction<T,Syn>& fraction<T,Syn>::operator/=(const fraction<T,Syn>& f) {
   auto_read_write_lock<const Syn> rwlock(f, *this);
   _numerator   *= f._denominator;
   _denominator *= f._numerator;
   return *this;
}

template <typename T, typename Syn>
template <typename U, typename S>
fraction<T,Syn>& fraction<T,Syn>::operator=(const fraction<U,S>& f) {
   auto_read_write_lock<const synchronizable> rwlock(f, *this);
   _numerator   = f._numerator;
   _denominator = f._denominator;
   return *this;
}

template <typename T, typename Syn>
template <typename U, typename S>
fraction<T,Syn>& fraction<T,Syn>::operator+=(const fraction<U,S>& f) {
   auto_read_write_lock<const synchronizable> rwlock(f, *this);
   _numerator = (_numerator * f._denominator) + (f._numerator * _denominator);
   _denominator *= f._denominator;
   return *this;
}

template <typename T, typename Syn>
template <typename U, typename S>
fraction<T,Syn>& fraction<T,Syn>::operator-=(const fraction<U,S>& f) {
   auto_read_write_lock<const synchronizable> rwlock(f, *this);
   _numerator = (_numerator * f._denominator) - (f._numerator * _denominator);
   _denominator *= f._denominator;
   return *this;
}

template <typename T, typename Syn>
template <typename U, typename S>
fraction<T,Syn>& fraction<T,Syn>::operator*=(const fraction<U,S>& f) {
   auto_read_write_lock<const synchronizable> rwlock(f, *this);
   _numerator   *= f._numerator;
   _denominator *= f._denominator;
   return *this;
}

template <typename T, typename Syn>
template <typename U, typename S>
fraction<T,Syn>& fraction<T,Syn>::operator/=(const fraction<U,S>& f) {
   auto_read_write_lock<const synchronizable> rwlock(f, *this);
   _numerator   *= f._denominator;
   _denominator *= f._numerator;
   return *this;
}

/***************************************************************************
 * Binary fraction-value operators.
 ***************************************************************************/

template <typename T, typename Syn>
inline fraction<T,Syn> operator+(const fraction<T,Syn>& f, const T& t) {
   return (fraction<T,Syn>(f) += t);
}

template <typename T, typename Syn>
inline fraction<T,Syn> operator-(const fraction<T,Syn>& f, const T& t) {
   return (fraction<T,Syn>(f) -= t);
}

template <typename T, typename Syn>
inline fraction<T,Syn> operator*(const fraction<T,Syn>& f, const T& t) {
   return (fraction<T,Syn>(f) *= t);
}

template <typename T, typename Syn>
inline fraction<T,Syn> operator/(const fraction<T,Syn>& f, const T& t) {
   return (fraction<T,Syn>(f) /= t);
}

template <typename T, typename Syn>
inline fraction<T,Syn> operator+(const T& t, const fraction<T,Syn>& f) {
   return (fraction<T,Syn>(f) += t);
}

template <typename T, typename Syn>
inline fraction<T,Syn> operator-(const T& t, const fraction<T,Syn>& f) {
   return (fraction<T,Syn>(t) -= f);
}

template <typename T, typename Syn>
inline fraction<T,Syn> operator*(const T& t, const fraction<T,Syn>& f) {
   return (fraction<T,Syn>(f) *= t);
}

template <typename T, typename Syn>
inline fraction<T,Syn> operator/(const T& t, const fraction<T,Syn>& f) {
   return (fraction<T,Syn>(t) /= f);
}

/***************************************************************************
 * Binary fraction-fraction operators.
 ***************************************************************************/

template <typename T, typename Syn>
inline fraction<T,Syn> operator+(
   const fraction<T,Syn>& f0, const fraction<T,Syn>& f1)
{
   return (fraction<T,Syn>(f0) += f1);
}

template <typename T, typename Syn>
inline fraction<T,Syn> operator-(
   const fraction<T,Syn>& f0, const fraction<T,Syn>& f1)
{
   return (fraction<T,Syn>(f0) -= f1);
}

template <typename T, typename Syn>
inline fraction<T,Syn> operator*(
   const fraction<T,Syn>& f0, const fraction<T,Syn>& f1)
{
   return (fraction<T,Syn>(f0) *= f1);
}

template <typename T, typename Syn>
inline fraction<T,Syn> operator/(
   const fraction<T,Syn>& f0, const fraction<T,Syn>& f1)
{
   return (fraction<T,Syn>(f0) /= f1);
}

/***************************************************************************
 * Binary fraction-value comparators.
 ***************************************************************************/

template <typename T, typename Syn>
bool operator==(const fraction<T,Syn>& f, const T& t) {
   auto_read_lock<const Syn> rlock(f);
   return (f._numerator == (t * f._denominator));
}

template <typename T, typename Syn>
bool operator!=(const fraction<T,Syn>& f, const T& t) {
   auto_read_lock<const Syn> rlock(f);
   return (f._numerator != (t * f._denominator));
}

template <typename T, typename Syn>
bool operator<(const fraction<T,Syn>& f, const T& t) {
   auto_read_lock<const Syn> rlock(f);
   return (
      (f._denominator < T()) ?
         (f._numerator > (t * f._denominator))
       : (f._numerator < (t * f._denominator))
   );
}

template <typename T, typename Syn>
bool operator>(const fraction<T,Syn>& f, const T& t) {
   auto_read_lock<const Syn> rlock(f);
   return (
      (f._denominator < T()) ?
         (f._numerator < (t * f._denominator))
       : (f._numerator > (t * f._denominator))
   );
}

template <typename T, typename Syn>
bool operator<=(const fraction<T,Syn>& f, const T& t) {
   auto_read_lock<const Syn> rlock(f);
   return (
      (f._denominator < T()) ?
         (f._numerator >= (t * f._denominator))
       : (f._numerator <= (t * f._denominator))
   );
}

template <typename T, typename Syn>
bool operator>=(const fraction<T,Syn>& f, const T& t) {
   auto_read_lock<const Syn> rlock(f);
   return (
      (f._denominator < T()) ?
         (f._numerator <= (t * f._denominator))
       : (f._numerator >= (t * f._denominator))
   );
}

template <typename T, typename Syn>
inline bool operator==(const T& t, const fraction<T,Syn>& f) {
   return (f == t);
}

template <typename T, typename Syn>
inline bool operator!=(const T& t, const fraction<T,Syn>& f) {
   return (f != t);
}

template <typename T, typename Syn>
inline bool operator<(const T& t, const fraction<T,Syn>& f) {
   return (f > t);
}

template <typename T, typename Syn>
inline bool operator>(const T& t, const fraction<T,Syn>& f) {
   return (f < t);
}

template <typename T, typename Syn>
inline bool operator<=(const T& t, const fraction<T,Syn>& f) {
   return (f >= t);
}

template <typename T, typename Syn>
inline bool operator>=(const T& t, const fraction<T,Syn>& f) {
   return (f <= t);
}

/***************************************************************************
 * Binary fraction-fraction comparators.
 ***************************************************************************/

template <typename T, typename Syn>
bool operator==(const fraction<T,Syn>& f0, const fraction<T,Syn>& f1) {
   auto_read_read_lock<const Syn> rrlock(f0, f1);
   return (
      (f0._numerator * f1._denominator) == (f1._numerator * f0._denominator)
   );
}

template <typename T, typename Syn>
bool operator!=(const fraction<T,Syn>& f0, const fraction<T,Syn>& f1) {
   auto_read_read_lock<const Syn> rrlock(f0, f1);
   return (
      (f0._numerator * f1._denominator) != (f1._numerator * f0._denominator)
   );
}

template <typename T, typename Syn>
bool operator<(const fraction<T,Syn>& f0, const fraction<T,Syn>& f1) {
   auto_read_read_lock<const Syn> rrlock(f0, f1);
   T n0d1 = f0._numerator * f1._denominator;
   T n1d0 = f1._numerator * f0._denominator;
   const T t_zero = T();
   bool d_neg =
      ((f0._denominator < t_zero) && (f1._denominator > t_zero)) ||
      ((f0._denominator > t_zero) && (f1._denominator < t_zero));
   return (d_neg ? (n0d1 > n1d0) : (n0d1 < n1d0));
}

template <typename T, typename Syn>
bool operator>(const fraction<T,Syn>& f0, const fraction<T,Syn>& f1) {
   auto_read_read_lock<const Syn> rrlock(f0, f1);
   T n0d1 = f0._numerator * f1._denominator;
   T n1d0 = f1._numerator * f0._denominator;
   const T t_zero = T();
   bool d_neg =
      ((f0._denominator < t_zero) && (f1._denominator > t_zero)) ||
      ((f0._denominator > t_zero) && (f1._denominator < t_zero));
   return (d_neg ? (n0d1 < n1d0) : (n0d1 > n1d0));
}

template <typename T, typename Syn>
bool operator<=(const fraction<T,Syn>& f0, const fraction<T,Syn>& f1) {
   auto_read_read_lock<const Syn> rrlock(f0, f1);
   T n0d1 = f0._numerator * f1._denominator;
   T n1d0 = f1._numerator * f0._denominator;
   const T t_zero = T();
   bool d_neg =
      ((f0._denominator < t_zero) && (f1._denominator > t_zero)) ||
      ((f0._denominator > t_zero) && (f1._denominator < t_zero));
   return (d_neg ? (n0d1 >= n1d0) : (n0d1 <= n1d0));
}

template <typename T, typename Syn>
bool operator>=(const fraction<T,Syn>& f0, const fraction<T,Syn>& f1) {
   auto_read_read_lock<const Syn> rrlock(f0, f1);
   T n0d1 = f0._numerator * f1._denominator;
   T n1d0 = f1._numerator * f0._denominator;
   const T t_zero = T();
   bool d_neg =
      ((f0._denominator < t_zero) && (f1._denominator > t_zero)) ||
      ((f0._denominator > t_zero) && (f1._denominator < t_zero));
   return (d_neg ? (n0d1 <= n1d0) : (n0d1 >= n1d0));
}

/***************************************************************************
 * Comparison method.
 ***************************************************************************/

/*
 * Comparison method.
 */
template <typename T, typename Syn>
int fraction<T,Syn>::compare_to(const fraction<T,Syn>& f) const {
   auto_read_read_lock<const Syn> rrlock(*this, f);
   T n0d1 = _numerator * f._denominator;
   T n1d0 = f._numerator * _denominator;
   const T t_zero = T();
   bool d_neg =
      ((_denominator < t_zero) && (f._denominator > t_zero)) ||
      ((_denominator > t_zero) && (f._denominator < t_zero));
   T t_diff = (d_neg ? (n1d0 - n0d1) : (n0d1 - n1d0));
   return ((t_diff > t_zero) ? 1 : ((t_diff < t_zero) ? -1 : 0));
}

/***************************************************************************
 * Unary operators.
 ***************************************************************************/

/*
 * Unary arithmetic operators.
 */
template <typename T, typename Syn>
fraction<T,Syn> fraction<T,Syn>::operator+() const {
   auto_read_lock<const Syn> rlock(*this);
   return fraction<T,Syn>(+(_numerator), _denominator);
}

template <typename T, typename Syn>
fraction<T,Syn> fraction<T,Syn>::operator-() const {
   auto_read_lock<const Syn> rlock(*this);
   return fraction<T,Syn>(-(_numerator), _denominator);
}

template <typename T, typename Syn>
fraction<T,Syn>& fraction<T,Syn>::operator++() {
   auto_write_lock<const Syn> wlock(*this);
   _numerator += _denominator;
   return *this;
}

template <typename T, typename Syn>
fraction<T,Syn>& fraction<T,Syn>::operator--() {
   auto_write_lock<const Syn> wlock(*this);
   _numerator -= _denominator;
   return *this;
}

template <typename T, typename Syn>
fraction<T,Syn> fraction<T,Syn>::operator++(int) {
   auto_write_lock<const Syn> wlock(*this);
   fraction<T,Syn> f(_numerator, _denominator);
   _numerator += _denominator;
   return f;
}

template <typename T, typename Syn>
fraction<T,Syn> fraction<T,Syn>::operator--(int) {
   auto_write_lock<const Syn> wlock(*this);
   fraction<T,Syn> f(_numerator, _denominator);
   _numerator -= _denominator;
   return f;
}

/*
 * Unary logical inversion operator.
 */
template <typename T, typename Syn>
bool fraction<T,Syn>::operator!() const {
   auto_read_lock<const Syn> rlock(*this);
   const T t_zero = T();
   return ((_numerator == t_zero) && (_denominator != t_zero));
}

/*
 * Unary multiplicative inversion operator.
 */
template <typename T, typename Syn>
fraction<T,Syn> fraction<T,Syn>::operator~() const {
   auto_read_lock<const Syn> rlock(*this);
   return fraction<T,Syn>(_denominator, _numerator);
}

/***************************************************************************
 * Sign and absolute value.
 ***************************************************************************/

/*
 * Return the sign (-1, 0, or 1) of the value represented by the fraction.
 */
template <typename T, typename Syn>
int fraction<T,Syn>::sign() const {
   auto_read_lock<const Syn> rlock(*this);
   const T t_zero = T();
   int sign_num =
      (_numerator > t_zero)   ? 1 : ((_numerator < t_zero)   ? -1 : 0);
   int sign_denom =
      (_denominator > t_zero) ? 1 : ((_denominator < t_zero) ? -1 : 0);
   return (sign_num * sign_denom);
}

/*
 * Absolute value.
 */
template <typename T, typename Syn>
fraction<T,Syn> abs(const fraction<T,Syn>& f) {
   auto_read_lock<const Syn> rlock(f);
   return fraction<T,Syn>(abs(f._numerator), abs(f._denominator));
}

/***************************************************************************
 * Numerator, denominator, and value.
 ***************************************************************************/

/*
 * Numerator.
 */
template <typename T, typename Syn>
T fraction<T,Syn>::numerator() const {
   auto_read_lock<const Syn> rlock(*this);
   return _numerator;
}

/*
 * Denominator.
 */
template <typename T, typename Syn>
T fraction<T,Syn>::denominator() const {
   auto_read_lock<const Syn> rlock(*this);
   return _denominator;
}

/*
 * Value.
 * Return the numerator divided by the denominator.
 */
template <typename T, typename Syn>
T fraction<T,Syn>::value() const {
   auto_read_lock<const Syn> rlock(*this);
   return (_numerator / _denominator);
}

/***************************************************************************
 * Reduction.
 ***************************************************************************/

/*
 * Reduce the fraction.
 * If T is an integer type, divide both the numerator and denominator by
 * their greatest common divisor.  Otherwise, if the denominator is nonzero,
 * divide the numerator by the denominator and set the denominator to one.
 */
template <typename T, typename Syn>
fraction<T,Syn>& fraction<T,Syn>::reduce() {
   auto_write_lock<const Syn> wlock(*this);
   /* check that denominator is nonzero */
   if (_denominator != T()) {
      /* compute gcd */
      const T t_one = T(1);
      T gcd = _numerator;
      T x = _denominator;
      do {
         T temp = x;
         x = gcd - ((gcd / x) * x);
         gcd = temp;
      } while (x >= t_one);
      /* divide */
      _numerator   /= gcd;
      _denominator /= gcd;
   }
   return *this;
}

} /* namespace math */

#endif
