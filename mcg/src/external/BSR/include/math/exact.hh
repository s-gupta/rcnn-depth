/*
 * Exact floating-point arithmetic (thread-safe).
 *
 * The exact class implements exact floating-point arithmetic (addition,
 * subtraction, and multiplication) using a multiple-term representation as 
 * described in:
 *
 * Jonathan Richard Shewchuk, "Adaptive Precision Floating-Point Arithmetic
 * and Fast Robust Geometric Predicates", Discrete & Computational Geometry
 * 18:305-363, 1997.
 *
 * Aribitrary precision significands are supported, however exponents are 
 * restricted to the usual range for double precision floating-point values.
 *
 * NOTE: Correct functioning of exact arithmetic requires that the processor's
 * internal floating-point registers not use extended precision bits.  Linking
 * a program with this module will cause all floating-point operations to 
 * round to double precision.
 */
#ifndef MATH__EXACT_HH
#define MATH__EXACT_HH

#include "concurrent/threads/synchronization/locks/auto_read_lock.hh"
#include "concurrent/threads/synchronization/locks/auto_read_read_lock.hh"
#include "concurrent/threads/synchronization/locks/auto_read_write_lock.hh"
#include "concurrent/threads/synchronization/locks/auto_write_lock.hh"
#include "concurrent/threads/synchronization/synchronizables/synchronizable.hh"
#include "concurrent/threads/synchronization/synchronizables/unsynchronized.hh"
#include "interfaces/comparable.hh"
#include "io/serialization/serial_input_stream.hh"
#include "io/serialization/serial_output_stream.hh"
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
using io::streams::ostream;
using lang::pointers::auto_ptr;

/*
 * Declare existence of exact_fp class.
 */
class exact_fp;

/*
 * Declare prototypes for exact_fp friend functions.
 */
ostream& operator<<(ostream&, const exact_fp&);

exact_fp operator+(const exact_fp&, const double&);
exact_fp operator-(const exact_fp&, const double&);
exact_fp operator*(const exact_fp&, const double&);

exact_fp operator+(const double&, const exact_fp&);
exact_fp operator-(const double&, const exact_fp&);
exact_fp operator*(const double&, const exact_fp&);

exact_fp operator+(const exact_fp&, const exact_fp&);
exact_fp operator-(const exact_fp&, const exact_fp&);
exact_fp operator*(const exact_fp&, const exact_fp&);

bool operator==(const exact_fp&, const double&);
bool operator!=(const exact_fp&, const double&);
bool operator< (const exact_fp&, const double&);
bool operator> (const exact_fp&, const double&);
bool operator<=(const exact_fp&, const double&);
bool operator>=(const exact_fp&, const double&);

bool operator==(const double&, const exact_fp&);
bool operator!=(const double&, const exact_fp&);
bool operator< (const double&, const exact_fp&);
bool operator> (const double&, const exact_fp&);
bool operator<=(const double&, const exact_fp&);
bool operator>=(const double&, const exact_fp&);

bool operator==(const exact_fp&, const exact_fp&);
bool operator!=(const exact_fp&, const exact_fp&);
bool operator< (const exact_fp&, const exact_fp&);
bool operator> (const exact_fp&, const exact_fp&);
bool operator<=(const exact_fp&, const exact_fp&);
bool operator>=(const exact_fp&, const exact_fp&);

exact_fp abs(const exact_fp&);

/*
 * Declare existence of exact class.
 */
template <typename Syn> class exact;

/*
 * Declare prototypes for exact template friend functions.
 */
template <typename Syn> ostream& operator<<(ostream&, const exact<Syn>&);

template <typename Syn> exact<Syn> operator+(const exact<Syn>&, const double&);
template <typename Syn> exact<Syn> operator-(const exact<Syn>&, const double&);
template <typename Syn> exact<Syn> operator*(const exact<Syn>&, const double&);

template <typename Syn> exact<Syn> operator+(const double&, const exact<Syn>&);
template <typename Syn> exact<Syn> operator-(const double&, const exact<Syn>&);
template <typename Syn> exact<Syn> operator*(const double&, const exact<Syn>&);

template <typename Syn>
exact<Syn> operator+(const exact<Syn>&, const exact<Syn>&);
template <typename Syn>
exact<Syn> operator-(const exact<Syn>&, const exact<Syn>&);
template <typename Syn>
exact<Syn> operator*(const exact<Syn>&, const exact<Syn>&);

template <typename Syn> bool operator==(const exact<Syn>&, const double&);
template <typename Syn> bool operator!=(const exact<Syn>&, const double&);
template <typename Syn> bool operator< (const exact<Syn>&, const double&);
template <typename Syn> bool operator> (const exact<Syn>&, const double&);
template <typename Syn> bool operator<=(const exact<Syn>&, const double&);
template <typename Syn> bool operator>=(const exact<Syn>&, const double&);

template <typename Syn> bool operator==(const double&, const exact<Syn>&);
template <typename Syn> bool operator!=(const double&, const exact<Syn>&);
template <typename Syn> bool operator< (const double&, const exact<Syn>&);
template <typename Syn> bool operator> (const double&, const exact<Syn>&);
template <typename Syn> bool operator<=(const double&, const exact<Syn>&);
template <typename Syn> bool operator>=(const double&, const exact<Syn>&);

template <typename Syn> bool operator==(const exact<Syn>&, const exact<Syn>&);
template <typename Syn> bool operator!=(const exact<Syn>&, const exact<Syn>&);
template <typename Syn> bool operator< (const exact<Syn>&, const exact<Syn>&);
template <typename Syn> bool operator> (const exact<Syn>&, const exact<Syn>&);
template <typename Syn> bool operator<=(const exact<Syn>&, const exact<Syn>&);
template <typename Syn> bool operator>=(const exact<Syn>&, const exact<Syn>&);

template <typename Syn> exact<Syn> abs(const exact<Syn>&);

/*
 * Exact floating-point base class.
 */
class exact_fp {
public:
   /*
    * Friend classes.
    */
   template <typename Syn> friend class exact;

   /*
    * Destructor.
    */
   ~exact_fp();
  
protected:
   /*
    * Constructors.
    */
   exact_fp();
   explicit exact_fp(const double&);

   /*
    * Copy constructor.
    */
   exact_fp(const exact_fp&);
   
   /*
    * Serialize.
    */
   void serialize(serial_output_stream&) const;

   /*
    * Deserialize.
    */
   static auto_ptr<exact_fp> deserialize(serial_input_stream&);
   
   /*
    * Formatted output to stream.
    */
   friend ostream& operator<<(ostream&, const exact_fp&);

   /*
    * Assignment operators: exact_fp-double.
    */
   exact_fp& operator=(const double&);
   exact_fp& operator+=(const double&);
   exact_fp& operator-=(const double&);
   exact_fp& operator*=(const double&);

   /*
    * Assignment operators: exact_fp-exact_fp.
    */
   exact_fp& operator=(const exact_fp&);
   exact_fp& operator+=(const exact_fp&);
   exact_fp& operator-=(const exact_fp&);
   exact_fp& operator*=(const exact_fp&);
   
   /*
    * Binary exact_fp-double operators.
    */
   friend exact_fp operator+(const exact_fp&, const double&);
   friend exact_fp operator-(const exact_fp&, const double&);
   friend exact_fp operator*(const exact_fp&, const double&);
   
   friend exact_fp operator+(const double&, const exact_fp&);
   friend exact_fp operator-(const double&, const exact_fp&);
   friend exact_fp operator*(const double&, const exact_fp&);

   /*
    * Binary exact-double operators.
    */
   template <typename Syn>
   friend exact<Syn> operator+(const exact<Syn>&, const double&);
   template <typename Syn>
   friend exact<Syn> operator-(const exact<Syn>&, const double&);
   template <typename Syn>
   friend exact<Syn> operator*(const exact<Syn>&, const double&);
   
   template <typename Syn>
   friend exact<Syn> operator+(const double&, const exact<Syn>&);
   template <typename Syn>
   friend exact<Syn> operator-(const double&, const exact<Syn>&);
   template <typename Syn>
   friend exact<Syn> operator*(const double&, const exact<Syn>&);

   /*
    * Binary exact_fp-exact_fp operators.
    */
   friend exact_fp operator+(const exact_fp&, const exact_fp&);
   friend exact_fp operator-(const exact_fp&, const exact_fp&);
   friend exact_fp operator*(const exact_fp&, const exact_fp&);
   
   /*
    * Binary exact-exact operators.
    */
   template <typename Syn>
   friend exact<Syn> operator+(const exact<Syn>&, const exact<Syn>&);
   template <typename Syn>
   friend exact<Syn> operator-(const exact<Syn>&, const exact<Syn>&);
   template <typename Syn>
   friend exact<Syn> operator*(const exact<Syn>&, const exact<Syn>&);

   /*
    * Binary exact_fp-double comparators.
    */
   friend bool operator==(const exact_fp&, const double&);
   friend bool operator!=(const exact_fp&, const double&);
   friend bool operator< (const exact_fp&, const double&);
   friend bool operator> (const exact_fp&, const double&);
   friend bool operator<=(const exact_fp&, const double&);
   friend bool operator>=(const exact_fp&, const double&);
   
   friend bool operator==(const double&, const exact_fp&);
   friend bool operator!=(const double&, const exact_fp&);
   friend bool operator< (const double&, const exact_fp&);
   friend bool operator> (const double&, const exact_fp&);
   friend bool operator<=(const double&, const exact_fp&);
   friend bool operator>=(const double&, const exact_fp&);

   /*
    * Binary exact-double comparators.
    */
   template <typename Syn>
   friend bool operator==(const exact<Syn>&, const double&);
   template <typename Syn>
   friend bool operator!=(const exact<Syn>&, const double&);
   template <typename Syn>
   friend bool operator< (const exact<Syn>&, const double&);
   template <typename Syn>
   friend bool operator> (const exact<Syn>&, const double&);
   template <typename Syn>
   friend bool operator<=(const exact<Syn>&, const double&);
   template <typename Syn>
   friend bool operator>=(const exact<Syn>&, const double&);
   
   template <typename Syn>
   friend bool operator==(const double&, const exact<Syn>&);
   template <typename Syn>
   friend bool operator!=(const double&, const exact<Syn>&);
   template <typename Syn>
   friend bool operator< (const double&, const exact<Syn>&);
   template <typename Syn>
   friend bool operator> (const double&, const exact<Syn>&);
   template <typename Syn>
   friend bool operator<=(const double&, const exact<Syn>&);
   template <typename Syn>
   friend bool operator>=(const double&, const exact<Syn>&);

   /*
    * Binary exact_fp-exact_fp comparators.
    */
   friend bool operator==(const exact_fp&, const exact_fp&);
   friend bool operator!=(const exact_fp&, const exact_fp&);
   friend bool operator< (const exact_fp&, const exact_fp&);
   friend bool operator> (const exact_fp&, const exact_fp&);
   friend bool operator<=(const exact_fp&, const exact_fp&);
   friend bool operator>=(const exact_fp&, const exact_fp&);

   /*
    * Binary exact-exact comparators.
    */
   template <typename Syn>
   friend bool operator==(const exact<Syn>&, const exact<Syn>&);
   template <typename Syn>
   friend bool operator!=(const exact<Syn>&, const exact<Syn>&);
   template <typename Syn>
   friend bool operator< (const exact<Syn>&, const exact<Syn>&);
   template <typename Syn>
   friend bool operator> (const exact<Syn>&, const exact<Syn>&);
   template <typename Syn>
   friend bool operator<=(const exact<Syn>&, const exact<Syn>&);
   template <typename Syn>
   friend bool operator>=(const exact<Syn>&, const exact<Syn>&);

   /*
    * Unary arithmetic operators.
    */
   exact_fp operator+() const;
   exact_fp operator-() const;
   exact_fp& operator++();
   exact_fp& operator--();
   exact_fp operator++(int);
   exact_fp operator--(int);

   /*
    * Unary logical inversion operator.
    */
   bool operator!() const;
   
   /*
    * Return the sign (-1, 0, or 1) of the exact value.
    */
   int sign() const;
   
   /*
    * Absolute value.
    */
   friend exact_fp abs(const exact_fp&);

   /*
    * Get size of representation (number of terms).
    */
   unsigned long size() const;

   /*
    * Compress the representation.
    * After compression, the most significant term is guaranteed to be a good
    * approximation to the true value.
    */
   exact_fp& compress();
   
   /*
    * Get floating-point approximation.
    * Return the approximation using the n most significant terms.
    * If n == 0, return the approximation using all terms.
    */
   double approx(unsigned long /* n */ = 0) const;
   
   /*
    * Return the relative error bound for double precision floating-point 
    * operations (2^-p where p is the number of bits in the significand).
    *
    *             |error(a O b)/(a O b)| <= epsilon
    *
    * where O denotes a floating-point operation (+, -, or *).
    */
   static double epsilon();

private:
   /************************************************************************
    * Data structures.
    ************************************************************************/
    
   unsigned long _size;       /* size of term list */
   unsigned long _size_alloc; /* allocated size */
   double* _fp;               /* term list (ordered by increasing magnitude) */

   /************************************************************************
    * Helper functions.
    ************************************************************************/
    
   /*
    * Swap the contents of two exact_fp objects.
    */
   static void swap(exact_fp&, exact_fp&);
   
   /*
    * Private constructor.
    * Allocate an uninitialized expansion of the given size.
    */
   exact_fp(unsigned long /* size */);

   /*
    * Split.
    * Given p-bit floating-point value x, return (floor(p/2))-bit value x_hi 
    * and (ceil(p/2)-1)-bit value x_low such that |x_hi| > |x_low| and 
    * x = x_hi + x_low is a nonoverlapping expansion.
    */
   static void split(double /* x */, double& /* x_hi */, double& /* x_low */);
   
   /*
    * Fast two-sum.
    * Given |a| > |b|, return nonoverlapping expansion x + y = a + b, where 
    * x is an approximation to a + b and y is the roundoff error.
    */
   static void fast_two_sum(
      double /* a */, double /* b */, double& /* x */, double& /* y */
   );
   
   /*
    * Fast two-diff.
    * Given |a| > |b|, return nonoverlapping expansion x + y = a - b, where 
    * x is an approximation to a - b and y is the roundoff error.
    */
   static void fast_two_diff(
      double /* a */, double /* b */, double& /* x */, double& /* y */
   );
  
   /*
    * Two-sum.
    * Given a, b, return nonoverlapping expansion x + y = a + b, where
    * x is an approximation to a + b and y is the roundoff error.
    */
   static void two_sum(
      double /* a */, double /* b */, double& /* x */, double& /* y */
   );
   
   /*
    * Two-diff.
    * Given a, b, return nonoverlapping expansion x + y = a - b, where
    * x is an approximation to a - b and y is the roundoff error.
    */
   static void two_diff(
      double /* a */, double /* b */, double& /* x */, double& /* y */
   );
   
   /*
    * Two-product.
    * Given a, b, return nonoverlapping expansion x + y = a * b, where
    * x is an approximation to a * b and y is the roundoff error.
    */
   static void two_product(
      double /* a */, double /* b */, double& /* x */, double& /* y */
   );
  
   /*
    * Grow expansion.
    * Given a nonoverlapping expansion and a value, return their sum as a 
    * nonoverlapping expansion.
    */
   static exact_fp grow_expansion(const exact_fp&, double);
   
   /*
    * Scale expansion.
    * Multiply the nonoverlapping expansion by the given value and return the
    * result as a nonoverlapping expansion.
    */
   static exact_fp scale_expansion(const exact_fp&, double);

   /*
    * Linear expansion sum.
    * Given two nonoverlapping expansions, return their sum as a nonoverlapping
    * expansion.
    */
   static exact_fp linear_expansion_sum(const exact_fp&, const exact_fp&);

   /*
    * Expansion product.
    * Given two nonoverlapping expansions, return their product as a
    * nonoverlapping expansion.
    */
   static exact_fp expansion_product(const exact_fp&, const exact_fp&);

   /*
    * Expansion product.
    * Given a nonoverlapping expansion e and a length N > 0 array of values v,
    * return the exact quantity:
    *
    *                            N-1
    *                            sum v[n]*e
    *                            n=0
    *
    * as a nonoverlapping expansion.
    */
   static exact_fp expansion_product(const exact_fp&, double*, unsigned long);
};

/*
 * Exact class.
 */
template <typename Syn = unsynchronized>
class exact : public comparable< exact<Syn> >,
              protected Syn {
public:
   /*
    * Friend classes.
    */
   template <typename S> friend class exact;

   /*
    * Constructors.
    */
   exact();
   explicit exact(const double&);

   /*
    * Copy constructors.
    */
   exact(const exact<Syn>&);
   template <typename S> exact(const exact<S>&);

   /*
    * Destructor.
    */
   virtual ~exact();

   /*
    * Serialize.
    */
   void serialize(serial_output_stream&) const;

   /*
    * Deserialize.
    */
   static auto_ptr< exact<Syn> > deserialize(serial_input_stream&);

   /*
    * Formatted output to stream.
    */
   friend ostream& operator<< <Syn>(ostream&, const exact<Syn>&);

   /*
    * Assignment operators: exact-double.
    */
   exact<Syn>& operator=(const double&);
   exact<Syn>& operator+=(const double&);
   exact<Syn>& operator-=(const double&);
   exact<Syn>& operator*=(const double&);

   /*
    * Assignment operators: exact-exact.
    */
   exact<Syn>& operator=(const exact<Syn>&);
   exact<Syn>& operator+=(const exact<Syn>&);
   exact<Syn>& operator-=(const exact<Syn>&);
   exact<Syn>& operator*=(const exact<Syn>&);

   template <typename S> exact<Syn>& operator=(const exact<S>&);
   template <typename S> exact<Syn>& operator+=(const exact<S>&);
   template <typename S> exact<Syn>& operator-=(const exact<S>&);
   template <typename S> exact<Syn>& operator*=(const exact<S>&);

   /*
    * Binary exact-double operators.
    */
   friend exact<Syn> operator+ <Syn>(const exact<Syn>&, const double&);
   friend exact<Syn> operator- <Syn>(const exact<Syn>&, const double&);
   friend exact<Syn> operator* <Syn>(const exact<Syn>&, const double&);
   
   friend exact<Syn> operator+ <Syn>(const double&, const exact<Syn>&);
   friend exact<Syn> operator- <Syn>(const double&, const exact<Syn>&);
   friend exact<Syn> operator* <Syn>(const double&, const exact<Syn>&);

   /*
    * Binary exact-exact operators.
    */
   friend exact<Syn> operator+ <Syn>(const exact<Syn>&, const exact<Syn>&);
   friend exact<Syn> operator- <Syn>(const exact<Syn>&, const exact<Syn>&);
   friend exact<Syn> operator* <Syn>(const exact<Syn>&, const exact<Syn>&);

   /*
    * Binary exact-double comparators.
    */
   friend bool operator== <Syn>(const exact<Syn>&, const double&);
   friend bool operator!= <Syn>(const exact<Syn>&, const double&);
   friend bool operator<  <Syn>(const exact<Syn>&, const double&);
   friend bool operator>  <Syn>(const exact<Syn>&, const double&);
   friend bool operator<= <Syn>(const exact<Syn>&, const double&);
   friend bool operator>= <Syn>(const exact<Syn>&, const double&);
   
   friend bool operator== <Syn>(const double&, const exact<Syn>&);
   friend bool operator!= <Syn>(const double&, const exact<Syn>&);
   friend bool operator<  <Syn>(const double&, const exact<Syn>&);
   friend bool operator>  <Syn>(const double&, const exact<Syn>&);
   friend bool operator<= <Syn>(const double&, const exact<Syn>&);
   friend bool operator>= <Syn>(const double&, const exact<Syn>&);

   /*
    * Binary exact-exact comparators.
    */
   friend bool operator== <Syn>(const exact<Syn>&, const exact<Syn>&);
   friend bool operator!= <Syn>(const exact<Syn>&, const exact<Syn>&);
   friend bool operator<  <Syn>(const exact<Syn>&, const exact<Syn>&);
   friend bool operator>  <Syn>(const exact<Syn>&, const exact<Syn>&);
   friend bool operator<= <Syn>(const exact<Syn>&, const exact<Syn>&);
   friend bool operator>= <Syn>(const exact<Syn>&, const exact<Syn>&);
   
   /*
    * Comparison method.
    */
   int compare_to(const exact<Syn>&) const;

   /*
    * Unary arithmetic operators.
    */
   exact<Syn> operator+() const;
   exact<Syn> operator-() const;
   exact<Syn>& operator++();
   exact<Syn>& operator--();
   exact<Syn> operator++(int);
   exact<Syn> operator--(int);

   /*
    * Unary logical inversion operator.
    */
   bool operator!() const;

   /*
    * Return the sign (-1, 0, or 1) of the exact value.
    */
   int sign() const;

   /*
    * Absolute value.
    */
   friend exact<Syn> abs<Syn>(const exact<Syn>&);

   /*
    * Get size of representation (number of terms).
    */
   unsigned long size() const;

   /*
    * Compress the representation.
    * After compression, the most significant term is guaranteed to be a good
    * approximation to the true value.
    */
   exact<Syn>& compress();

   /*
    * Get floating-point approximation.
    * Return the approximation using the n most significant terms.
    * If n == 0, return the approximation using all terms.
    */
   double approx(unsigned long /* n */ = 0) const;

   /*
    * Return the relative error bound for double precision floating-point 
    * operations (2^-p where p is the number of bits in the significand).
    *
    *             |error(a O b)/(a O b)| <= epsilon
    *
    * where O denotes a floating-point operation (+, -, or *).
    */
   static double epsilon();

protected:
   /*
    * Convert an exact_fp object into an exact<Syn> object.
    * The original exact_fp object is destroyed in the process.
    */
   static exact<Syn> to_exact(exact_fp&);

private:
   exact_fp _e_fp;   /* underlying exact_fp object */
};

/***************************************************************************
 * Constructors and destructor (exact).
 ***************************************************************************/

/*
 * Default constructor.
 * Return the exact value zero.
 */
template <typename Syn>
exact<Syn>::exact()
 : Syn(),
   _e_fp()
{ }

/*
 * Constructor.
 * Return the given value.
 */
template <typename Syn>
exact<Syn>::exact(const double& x)
 : Syn(),
   _e_fp(x)
{ }

/*
 * Copy constructor.
 */
template <typename Syn>
exact<Syn>::exact(const exact<Syn>& e)
 : Syn(),
   _e_fp()
{
   auto_read_lock<const Syn> rlock(e);
   exact_fp e_copy(e._e_fp);
   exact_fp::swap(_e_fp, e_copy);
}

/*
 * Type conversion copy constructor.
 */
template <typename Syn>
template <typename S>
exact<Syn>::exact(const exact<S>& e)
 : Syn(),
   _e_fp()
{
   auto_read_lock<const S> rlock(e);
   exact_fp e_copy(e._e_fp);
   exact_fp::swap(_e_fp, e_copy);
}

/*
 * Destructor.
 */
template <typename Syn>
exact<Syn>::~exact() {
   /* do nothing */
}

/***************************************************************************
 * Serialization (exact).
 ***************************************************************************/

/*
 * Serialize.
 */
template <typename Syn>
void exact<Syn>::serialize(serial_output_stream& s) const {
   auto_read_lock<const Syn> rlock(*this);
   _e_fp.serialize(s);
}

/*
 * Deserialize.
 */
template <typename Syn>
auto_ptr< exact<Syn> > exact<Syn>::deserialize(serial_input_stream& s) {
   auto_ptr<exact_fp> e_fp = exact_fp::deserialize(s);
   auto_ptr< exact<Syn> > e(new exact<Syn>());
   e._e_fp = *e_fp;
   return e;
}

/***************************************************************************
 * I/O (exact).
 ***************************************************************************/

/*
 * Formatted output to stream.
 */
template <typename Syn>
ostream& operator<<(ostream& os, const exact<Syn>& e) {
   auto_read_lock<const Syn> rlock(e);
   return (os << e._e_fp);
}

/***************************************************************************
 * Assignment operators: exact-double.
 ***************************************************************************/

template <typename Syn>
exact<Syn>& exact<Syn>::operator=(const double& x) {
   auto_write_lock<const Syn> wlock(*this);
   _e_fp = x;
   return *this;
}

template <typename Syn>
exact<Syn>& exact<Syn>::operator+=(const double& x) {
   auto_write_lock<const Syn> wlock(*this);
   _e_fp += x;
   return *this;
}

template <typename Syn>
exact<Syn>& exact<Syn>::operator-=(const double& x) {
   auto_write_lock<const Syn> wlock(*this);
   _e_fp -= x;
   return *this;
}

template <typename Syn>
exact<Syn>& exact<Syn>::operator*=(const double& x) {
   auto_write_lock<const Syn> wlock(*this);
   _e_fp *= x;
   return *this;
}

/***************************************************************************
 * Assignment operators: exact-exact.
 ***************************************************************************/

template <typename Syn>
exact<Syn>& exact<Syn>::operator=(const exact<Syn>& e) {
   auto_read_write_lock<const Syn> rwlock(e, *this);
   _e_fp = e._e_fp;
   return *this;
}

template <typename Syn>
exact<Syn>& exact<Syn>::operator+=(const exact<Syn>& e) {
   auto_read_write_lock<const Syn> rwlock(e, *this);
   _e_fp += e._e_fp;
   return *this;
}

template <typename Syn>
exact<Syn>& exact<Syn>::operator-=(const exact<Syn>& e) {
   auto_read_write_lock<const Syn> rwlock(e, *this);
   _e_fp -= e._e_fp;
   return *this;
}

template <typename Syn>
exact<Syn>& exact<Syn>::operator*=(const exact<Syn>& e) {
   auto_read_write_lock<const Syn> rwlock(e, *this);
   _e_fp *= e._e_fp;
   return *this;
}

template <typename Syn>
template <typename S>
exact<Syn>& exact<Syn>::operator=(const exact<S>& e) {
   auto_read_write_lock<const synchronizable> rwlock(e, *this);
   _e_fp = e._e_fp;
   return *this;
}

template <typename Syn>
template <typename S>
exact<Syn>& exact<Syn>::operator+=(const exact<S>& e) {
   auto_read_write_lock<const synchronizable> rwlock(e, *this);
   _e_fp += e._e_fp;
   return *this;
}

template <typename Syn>
template <typename S>
exact<Syn>& exact<Syn>::operator-=(const exact<S>& e) {
   auto_read_write_lock<const synchronizable> rwlock(e, *this);
   _e_fp -= e._e_fp;
   return *this;
}

template <typename Syn>
template <typename S>
exact<Syn>& exact<Syn>::operator*=(const exact<S>& e) {
   auto_read_write_lock<const synchronizable> rwlock(e, *this);
   _e_fp *= e._e_fp;
   return *this;
}

/***************************************************************************
 * Binary exact-double operators.
 ***************************************************************************/

template <typename Syn>
exact<Syn> operator+(const exact<Syn>& e, const double& x) {
   auto_read_lock<const Syn> rlock(e);
   exact_fp e_result = e._e_fp + x;
   return exact<Syn>::to_exact(e_result);
}

template <typename Syn>
exact<Syn> operator-(const exact<Syn>& e, const double& x) {
   auto_read_lock<const Syn> rlock(e);
   exact_fp e_result = e._e_fp - x;
   return exact<Syn>::to_exact(e_result);
}

template <typename Syn>
exact<Syn> operator*(const exact<Syn>& e, const double& x) {
   auto_read_lock<const Syn> rlock(e);
   exact_fp e_result = e._e_fp * x;
   return exact<Syn>::to_exact(e_result);
}

template <typename Syn>
exact<Syn> operator+(const double& x, const exact<Syn>& e) {
   auto_read_lock<const Syn> rlock(e);
   exact_fp e_result = x + e._e_fp;
   return exact<Syn>::to_exact(e_result);
}

template <typename Syn>
exact<Syn> operator-(const double& x, const exact<Syn>& e) {
   auto_read_lock<const Syn> rlock(e);
   exact_fp e_result = x - e._e_fp;
   return exact<Syn>::to_exact(e_result);
}

template <typename Syn>
exact<Syn> operator*(const double& x, const exact<Syn>& e) {
   auto_read_lock<const Syn> rlock(e);
   exact_fp e_result = x * e._e_fp;
   return exact<Syn>::to_exact(e_result);
}

/***************************************************************************
 * Binary exact-exact operators.
 ***************************************************************************/

template <typename Syn>
exact<Syn> operator+(const exact<Syn>& e0, const exact<Syn>& e1) {
   auto_read_read_lock<const Syn> rrlock(e0, e1);
   exact_fp e = e0._e_fp + e1._e_fp;
   return exact<Syn>::to_exact(e);
}

template <typename Syn>
exact<Syn> operator-(const exact<Syn>& e0, const exact<Syn>& e1) {
   auto_read_read_lock<const Syn> rrlock(e0, e1);
   exact_fp e = e0._e_fp - e1._e_fp;
   return exact<Syn>::to_exact(e);
}

template <typename Syn>
exact<Syn> operator*(const exact<Syn>& e0, const exact<Syn>& e1) {
   auto_read_read_lock<const Syn> rrlock(e0, e1);
   exact_fp e = e0._e_fp * e1._e_fp;
   return exact<Syn>::to_exact(e);
}

/***************************************************************************
 * Binary exact-double comparators.
 ***************************************************************************/

template <typename Syn>
bool operator==(const exact<Syn>& e, const double& x) {
   auto_read_lock<const Syn> rlock(e);
   return (e._e_fp == x);
}

template <typename Syn>
bool operator!=(const exact<Syn>& e, const double& x) {
   auto_read_lock<const Syn> rlock(e);
   return (e._e_fp != x);
}

template <typename Syn>
bool operator<(const exact<Syn>& e, const double& x) {
   auto_read_lock<const Syn> rlock(e);
   return (e._e_fp < x);
}

template <typename Syn>
bool operator>(const exact<Syn>& e, const double& x) {
   auto_read_lock<const Syn> rlock(e);
   return (e._e_fp > x);
}

template <typename Syn>
bool operator<=(const exact<Syn>& e, const double& x) {
   auto_read_lock<const Syn> rlock(e);
   return (e._e_fp <= x);
}

template <typename Syn>
bool operator>=(const exact<Syn>& e, const double& x) {
   auto_read_lock<const Syn> rlock(e);
   return (e._e_fp >= x);
}

template <typename Syn>
bool operator==(const double& x, const exact<Syn>& e) {
   auto_read_lock<const Syn> rlock(e);
   return (x == e._e_fp);
}

template <typename Syn>
bool operator!=(const double& x, const exact<Syn>& e) {
   auto_read_lock<const Syn> rlock(e);
   return (x != e._e_fp);
}

template <typename Syn>
bool operator<(const double& x, const exact<Syn>& e) {
   auto_read_lock<const Syn> rlock(e);
   return (x < e._e_fp);
}

template <typename Syn>
bool operator>(const double& x, const exact<Syn>& e) {
   auto_read_lock<const Syn> rlock(e);
   return (x > e._e_fp);
}

template <typename Syn>
bool operator<=(const double& x, const exact<Syn>& e) {
   auto_read_lock<const Syn> rlock(e);
   return (x <= e._e_fp);
}

template <typename Syn>
bool operator>=(const double& x, const exact<Syn>& e) {
   auto_read_lock<const Syn> rlock(e);
   return (x >= e._e_fp);
}

/***************************************************************************
 * Binary exact-exact comparators.
 ***************************************************************************/

template <typename Syn>
bool operator==(const exact<Syn>& e0, const exact<Syn>& e1) {
   auto_read_read_lock<const Syn> rrlock(e0, e1);
   return (e0._e_fp == e1._e_fp);
}

template <typename Syn>
bool operator!=(const exact<Syn>& e0, const exact<Syn>& e1) {
   auto_read_read_lock<const Syn> rrlock(e0, e1);
   return (e0._e_fp != e1._e_fp);
}

template <typename Syn>
bool operator<(const exact<Syn>& e0, const exact<Syn>& e1) {
   auto_read_read_lock<const Syn> rrlock(e0, e1);
   return (e0._e_fp < e1._e_fp);
}

template <typename Syn>
bool operator>(const exact<Syn>& e0, const exact<Syn>& e1) {
   auto_read_read_lock<const Syn> rrlock(e0, e1);
   return (e0._e_fp > e1._e_fp);
}

template <typename Syn>
bool operator<=(const exact<Syn>& e0, const exact<Syn>& e1) {
   auto_read_read_lock<const Syn> rrlock(e0, e1);
   return (e0._e_fp <= e1._e_fp);
}

template <typename Syn>
bool operator>=(const exact<Syn>& e0, const exact<Syn>& e1) {
   auto_read_read_lock<const Syn> rrlock(e0, e1);
   return (e0._e_fp >= e1._e_fp);
}

/***************************************************************************
 * Comparison method.
 ***************************************************************************/

/*
 * Comparison method.
 */
template <typename Syn>
int exact<Syn>::compare_to(const exact<Syn>& e) const {
   auto_read_read_lock<const Syn> rrlock(*this, e);
   exact_fp e_diff = _e_fp - e._e_fp;
   return e_diff.sign();
}

/***************************************************************************
 * Unary operators (exact).
 ***************************************************************************/

/*
 * Unary arithmetic operators.
 */
template <typename Syn>
exact<Syn> exact<Syn>::operator+() const {
   auto_read_lock<const Syn> rlock(*this);
   exact_fp e = +_e_fp;
   return exact<Syn>::to_exact(e);
}

template <typename Syn>
exact<Syn> exact<Syn>::operator-() const {
   auto_read_lock<const Syn> rlock(*this);
   exact_fp e = -_e_fp;
   return exact<Syn>::to_exact(e);
}

template <typename Syn>
exact<Syn>& exact<Syn>::operator++() {
   auto_write_lock<const Syn> wlock(*this);
   ++(_e_fp);
   return *this;
}

template <typename Syn>
exact<Syn>& exact<Syn>::operator--() {
   auto_write_lock<const Syn> wlock(*this);
   --(_e_fp);
   return *this;
}

template <typename Syn>
exact<Syn> exact<Syn>::operator++(int) {
   auto_write_lock<const Syn> wlock(*this);
   exact_fp e = _e_fp++;
   return exact<Syn>::to_exact(e);
}

template <typename Syn>
exact<Syn> exact<Syn>::operator--(int) {
   auto_write_lock<const Syn> wlock(*this);
   exact_fp e = _e_fp--;
   return exact<Syn>::to_exact(e);
}

/*
 * Unary logical inversion operator.
 */
template <typename Syn>
bool exact<Syn>::operator!() const {
   auto_read_lock<const Syn> rlock(*this);
   return !_e_fp;
}

/***************************************************************************
 * Sign and absolute value (exact).
 ***************************************************************************/

/*
 * Return the sign (-1, 0, or 1) of the exact value.
 */
template <typename Syn>
int exact<Syn>::sign() const {
   auto_read_lock<const Syn> rlock(*this);
   return _e_fp.sign();
}

/*
 * Absolute value.
 */
template <typename Syn>
exact<Syn> abs(const exact<Syn>& e) {
   auto_read_lock<const Syn> rlock(e);
   exact_fp e_abs = abs(e._e_fp);
   return exact<Syn>::to_exact(e_abs);
}

/***************************************************************************
 * Size and approximation (exact).
 ***************************************************************************/

/*
 * Get size of representation (number of terms).
 */
template <typename Syn>
unsigned long exact<Syn>::size() const {
   auto_read_lock<const Syn> rlock(*this);
   return _e_fp.size();
}

/*
 * Compress the representation.
 * After compression, the most significant term is guaranteed to be a good
 * approximation to the true value.
 */
template <typename Syn>
exact<Syn>& exact<Syn>::compress() {
   auto_write_lock<const Syn> wlock(*this);
   _e_fp.compress();
   return *this;
}

/*
 * Get floating-point approximation.
 * Return the approximation using the n most significant terms.
 * If n == 0, return the approximation using all terms.
 */
template <typename Syn>
double exact<Syn>::approx(unsigned long n_terms) const {
   auto_read_lock<const Syn> rlock(*this);
   return _e_fp.approx(n_terms);
}

/*
 * Return the relative error bound for double precision floating-point 
 * operations (2^-p where p is the number of bits in the significand).
 *
 *             |error(a O b)/(a O b)| <= epsilon
 *
 * where O denotes a floating-point operation (+, -, or *).
 */
template <typename Syn>
double exact<Syn>::epsilon() {
   return exact_fp::epsilon();
}

/***************************************************************************
 * Conversion to exact.
 ***************************************************************************/

/*
 * Convert an exact_fp object into an exact<Syn> object.
 * The original exact_fp object is destroyed in the process.
 */
template <typename Syn>
exact<Syn> exact<Syn>::to_exact(exact_fp& e) {
   exact<Syn> e_syn;
   exact_fp::swap(e_syn._e_fp, e);
   return e_syn;
}

} /* namespace math */

#endif
