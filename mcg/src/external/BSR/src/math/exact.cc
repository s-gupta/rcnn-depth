/*
 * Exact floating-point arithmetic.
 */
#include "io/serialization/serial_input_stream.hh"
#include "io/serialization/serial_output_stream.hh"
#include "io/streams/ostream.hh"
#include "lang/pointers/auto_ptr.hh"
#include "math/exact.hh"
#include "math/math.hh"

namespace math {
/*
 * Imports.
 */
using io::serialization::serial_input_stream;
using io::serialization::serial_output_stream;
using io::streams::ostream;
using lang::pointers::auto_ptr;

/***************************************************************************
 * Initialization and error bounds.
 ***************************************************************************/

namespace {
   /*
    * Initialize the floating-point unit for exact arithmetic by setting it 
    * to use exactly double precision (no extended precision bits).
    *
    * Also return the relative error bound and splitter for exact double
    * precision floating-point arithmetic.
    */
   void exact_fp_init(double& epsilon, double& splitter) {
      /* compute relative error bound and splitter */
      epsilon  = 1.0;
      splitter = 1.0;
      double approx      = 1.0;
      double approx_prev = 1.0;
      bool increment_splitter = true;
      do {
         approx_prev = approx;
         epsilon = 0.5 * epsilon;
         approx  = 1.0 + epsilon;
         if (increment_splitter)
            splitter *= 2.0;
         increment_splitter = !increment_splitter;
      } while ((approx != 1.0) && (approx != approx_prev));
      splitter += 1.0;
   }

   /*
    * Compute and return the relative error bound for floating-point
    * operations (2^-p where p is the number of bits in the significand).
    */
   double exact_fp_compute_epsilon() {
      double epsilon  = 0;
      double splitter = 0;
      exact_fp_init(epsilon, splitter);
      return epsilon;
   }
   
   /*
    * Compute and return 2^(ceiling(p/2)) + 1, where p is the number of bits
    * in a double precision floating-point significand.
    */
   double exact_fp_compute_splitter() {
      double epsilon  = 0;
      double splitter = 0;
      exact_fp_init(epsilon, splitter);
      return splitter;
   }
   
   /*
    * Global static variables for storing the relative error bound and
    * splitting value.
    *
    * Note that global variables are required here (instead of local statics
    * set with the initialize on first use idiom) to guarantee that the 
    * floating-point unit is initialized before everything else.
    *
    * All other static program variables which may interact with the exact
    * class should use the initialize on first use idiom.
    */
   static const double exact_fp_epsilon  = exact_fp_compute_epsilon();

   static const double exact_fp_splitter = exact_fp_compute_splitter();
} /* namespace */

/*
 * Return the relative error bound for double precision floating-point 
 * operations (2^-p where p is the number of bits in the significand).
 *
 *             error(a O b)/(a O b) <= epsilon
 *
 * where O denotes a floating-point operation (+, -, or *).
 */
double exact_fp::epsilon() {
   return exact_fp_epsilon;
}

/***************************************************************************
 * Constructors and destructor (exact_fp).
 ***************************************************************************/

/*
 * Default constructor.
 * Return the exact value zero.
 */
exact_fp::exact_fp()
 : _size(0),
   _size_alloc(0),
   _fp(NULL)
{ }

/*
 * Constructor.
 * Return the given value.
 */
exact_fp::exact_fp(const double& x)
 : _size(1),
   _size_alloc(1),
   _fp(new double[1]())
{
   _fp[0] = x;
   if (x == 0)
      _size = 0;
}

/*
 * Private constructor.
 * Allocate an uninitialized expansion of the given size.
 */
exact_fp::exact_fp(unsigned long size)
 : _size(size),
   _size_alloc(size),
   _fp(new double[size]())
{ }

/*
 * Copy constructor.
 */
exact_fp::exact_fp(const exact_fp& e)
 : _size(e._size),
   _size_alloc(_size),
   _fp((_size > 0) ? (new double[_size]()) : NULL)
{
   for (unsigned long n = 0; n < _size; n++)
      _fp[n] = e._fp[n];
}

/*
 * Destructor.
 */
exact_fp::~exact_fp() {
   delete [] _fp;
}

/***************************************************************************
 * Serialization (exact_fp).
 ***************************************************************************/

/*
 * Serialize.
 */
void exact_fp::serialize(serial_output_stream& s) const {
   s << _size;
   for (unsigned long n = 0; n < _size; n++)
      s << _fp[n];
}

/*
 * Deserialize.
 */
auto_ptr<exact_fp> exact_fp::deserialize(serial_input_stream& s) {
   unsigned long size;
   s >> size;
   auto_ptr<exact_fp> e(new exact_fp(size));
   for (unsigned long n = 0; n < size; n++)
      s >> e->_fp[n];
   return e;
}

/***************************************************************************
 * I/O (exact_fp).
 ***************************************************************************/

/*
 * Formatted output to stream.
 */
ostream& operator<<(ostream& os, const exact_fp& e) {
   os << e.approx();
   return os;
}

/***************************************************************************
 * Assignment operators: exact_fp-double.
 ***************************************************************************/

exact_fp& exact_fp::operator=(const double& x) {
   if (x == 0) {
      _size = 0;
   } else {
      if (_size_alloc == 0) {
         _fp = new double[1]();
         _size_alloc = 1;
      }
      _size = 1;
      _fp[0] = x;
   }
   return *this;
}

exact_fp& exact_fp::operator+=(const double& x) {
   exact_fp e_sum = exact_fp::grow_expansion(*this, x);
   exact_fp::swap(*this, e_sum);
   return *this;
}

exact_fp& exact_fp::operator-=(const double& x) {
   exact_fp e_diff = exact_fp::grow_expansion(*this, -x);
   exact_fp::swap(*this, e_diff);
   return *this;
}

exact_fp& exact_fp::operator*=(const double& x) {
   exact_fp e_prod = exact_fp::scale_expansion(*this, x);
   exact_fp::swap(*this, e_prod);
   return *this;
}

/***************************************************************************
 * Assignment operators: exact_fp-exact_fp.
 ***************************************************************************/

exact_fp& exact_fp::operator=(const exact_fp& e) {
   if (_size_alloc < e._size) {
      double* fp = new double[e._size]();
      delete _fp;
      _fp = fp;
      _size_alloc = e._size;
   }
   _size = e._size;
   for (unsigned long n = 0; n < _size; n++)
      _fp[n] = e._fp[n];
   return *this;
}

exact_fp& exact_fp::operator+=(const exact_fp& e) {
   exact_fp e_sum = exact_fp::linear_expansion_sum(*this, e);
   exact_fp::swap(*this, e_sum);
   return *this;
}

exact_fp& exact_fp::operator-=(const exact_fp& e) {
   exact_fp e_diff = exact_fp::linear_expansion_sum(*this, -e);
   exact_fp::swap(*this, e_diff);
   return *this;
}

exact_fp& exact_fp::operator*=(const exact_fp& e) {
   exact_fp e_prod = exact_fp::expansion_product(*this, e);
   exact_fp::swap(*this, e_prod);
   return *this;
}

/***************************************************************************
 * Binary exact_fp-double operators.
 ***************************************************************************/

exact_fp operator+(const exact_fp& e, const double& x) {
   return exact_fp::grow_expansion(e, x);
}

exact_fp operator-(const exact_fp& e, const double& x) {
   return exact_fp::grow_expansion(e, -x);
}

exact_fp operator*(const exact_fp& e, const double& x) {
   return exact_fp::scale_expansion(e, x);
}

exact_fp operator+(const double& x, const exact_fp& e) {
   return exact_fp::grow_expansion(e, x);
}

exact_fp operator-(const double& x, const exact_fp& e) {
   exact_fp e_diff = exact_fp::grow_expansion(e, -x);
   for (unsigned long n = 0; n < e_diff._size; n++)
      e_diff._fp[n] = -(e_diff._fp[n]);
   return e_diff;
}

exact_fp operator*(const double& x, const exact_fp& e) {
   return exact_fp::scale_expansion(e, x);
}

/***************************************************************************
 * Binary exact_fp-exact_fp operators.
 ***************************************************************************/

exact_fp operator+(const exact_fp& e0, const exact_fp& e1) {
   return exact_fp::linear_expansion_sum(e0, e1);
}

exact_fp operator-(const exact_fp& e0, const exact_fp& e1) {
   return exact_fp::linear_expansion_sum(e0, -e1);
}

exact_fp operator*(const exact_fp& e0, const exact_fp& e1) {
   return exact_fp::expansion_product(e0, e1);
}

/***************************************************************************
 * Binary exact_fp-double comparators.
 ***************************************************************************/

bool operator==(const exact_fp& e, const double& x) {
   exact_fp e_diff = e - x;
   return (e_diff.sign() == 0);
}

bool operator!=(const exact_fp& e, const double& x) {
   exact_fp e_diff = e - x;
   return (e_diff.sign() != 0);
}

bool operator<(const exact_fp& e, const double& x) {
   exact_fp e_diff = e - x;
   return (e_diff.sign() < 0);
}

bool operator>(const exact_fp& e, const double& x) {
   exact_fp e_diff = e - x;
   return (e_diff.sign() > 0);
}

bool operator<=(const exact_fp& e, const double& x) {
   exact_fp e_diff = e - x;
   return (e_diff.sign() <= 0);
}

bool operator>=(const exact_fp& e, const double& x) {
   exact_fp e_diff = e - x;
   return (e_diff.sign() >= 0);
}

bool operator==(const double& x, const exact_fp& e) {
   return (e == x);
}

bool operator!=(const double& x, const exact_fp& e) {
   return (e != x);
}

bool operator<(const double& x, const exact_fp& e) {
   return (e > x);
}

bool operator>(const double& x, const exact_fp& e) {
   return (e < x);
}

bool operator<=(const double& x, const exact_fp& e) {
   return (e >= x);
}

bool operator>=(const double& x, const exact_fp& e) {
   return (e <= x);
}

/***************************************************************************
 * Binary exact_fp-exact_fp comparators.
 ***************************************************************************/

bool operator==(const exact_fp& e0, const exact_fp& e1) {
   exact_fp e_diff = e0 - e1;
   return (e_diff.sign() == 0);
}

bool operator!=(const exact_fp& e0, const exact_fp& e1) {
   exact_fp e_diff = e0 - e1;
   return (e_diff.sign() != 0);
}

bool operator<(const exact_fp& e0, const exact_fp& e1) {
   exact_fp e_diff = e0 - e1;
   return (e_diff.sign() < 0);
}

bool operator>(const exact_fp& e0, const exact_fp& e1) {
   exact_fp e_diff = e0 - e1;
   return (e_diff.sign() > 0);
}

bool operator<=(const exact_fp& e0, const exact_fp& e1) {
   exact_fp e_diff = e0 - e1;
   return (e_diff.sign() <= 0);
}

bool operator>=(const exact_fp& e0, const exact_fp& e1) {
   exact_fp e_diff = e0 - e1;
   return (e_diff.sign() >= 0);
}

/***************************************************************************
 * Unary operators (exact_fp).
 ***************************************************************************/

/*
 * Unary arithmetic operators.
 */
exact_fp exact_fp::operator+() const {
   return *this;
}

exact_fp exact_fp::operator-() const {
   exact_fp e(_size);
   for (unsigned long n = 0; n < _size; n++)
      e._fp[n] = -(_fp[n]);
   return e;
}

exact_fp& exact_fp::operator++() {
   return (*this += 1.0);
}

exact_fp& exact_fp::operator--() {
   return (*this -= 1.0);
}

exact_fp exact_fp::operator++(int) {
   exact_fp e(*this);
   *this += 1.0;
   return e;
}

exact_fp exact_fp::operator--(int) {
   exact_fp e(*this);
   *this -= 1.0;
   return e;
}

/*
 * Unary logical inversion operator.
 */
bool exact_fp::operator!() const {
   return (_size == 0);
}

/***************************************************************************
 * Sign and absolute value (exact_fp).
 ***************************************************************************/

/*
 * Return the sign (-1, 0, or 1) of the exact value.
 */
int exact_fp::sign() const {
   if (_size == 0)
      return 0;
   else if (_fp[_size - 1] > 0)
      return 1;
   else
      return -1;
}

/*
 * Absolute value.
 */
exact_fp abs(const exact_fp& e) {
   return ((e.sign() < 0) ? (-e) : e);
}

/***************************************************************************
 * Size and approximation (exact_fp).
 ***************************************************************************/

/*
 * Get size of representation (number of terms).
 */
unsigned long exact_fp::size() const {
   return _size;
}

/*
 * Compress the representation.
 * After compression, the most significant term is guaranteed to be a good
 * approximation to the true value.
 */
exact_fp& exact_fp::compress() {
   /* check if expansion is trivial */
   if (_size <= 1)
      return *this;
   /* traverse from largest to smallest component, replacing adjacent pairs */
   double Q = _fp[_size - 1];
   double Q_next = 0, q = 0;
   unsigned long bottom = _size - 1;
   for (unsigned long n = 1; n < _size; n++) {
      unsigned long pos = _size - 1 - n;
      exact_fp::fast_two_sum(Q, _fp[pos], Q_next, q);
      if (q != 0) {
         _fp[bottom--] = Q_next;
         Q = q;
      } else {
         Q = Q_next;
      }
   }
   _fp[bottom] = Q;
   /* traverse from smallest to largest component, clipping overlapping bits */
   unsigned long top = 0;
   for (unsigned long n = bottom + 1; n < _size; n++) {
      exact_fp::fast_two_sum(_fp[n], Q, Q_next, q);
      Q = Q_next;
      if (q != 0)
         _fp[top++] = q;
   }
   _fp[top++] = Q;
   _size = top;
   return *this;
}

/*
 * Get floating-point approximation.
 * Return the approximation using the n most significant terms.
 * If n == 0, return the approximation using all terms.
 */
double exact_fp::approx(unsigned long n_terms) const {
   /* check number of terms */
   if ((n_terms == 0) || (n_terms > _size))
      n_terms = _size;
   /* sum terms */
   double x = 0;
   for (unsigned long n = (_size - n_terms); n < _size; n++)
      x += _fp[n];
   return x;
}

/***************************************************************************
 * Helper functions (exact_fp).
 ***************************************************************************/

/*
 * Swap the contents of two exact_fp objects.
 */
void exact_fp::swap(exact_fp& e0, exact_fp& e1) {
   unsigned long temp_size       = e0._size;
   unsigned long temp_size_alloc = e0._size_alloc;
   double*       temp_fp         = e0._fp;
   e0._size       = e1._size;
   e0._size_alloc = e1._size_alloc;
   e0._fp         = e1._fp;
   e1._size       = temp_size;
   e1._size_alloc = temp_size_alloc;
   e1._fp         = temp_fp;
}

/*
 * Split.
 * Given p-bit floating-point value x, return (floor(p/2))-bit value x_hi 
 * and (ceil(p/2)-1)-bit value x_low such that |x_hi| > |x_low| and 
 * x = x_hi + x_low is a nonoverlapping expansion.
 */
void exact_fp::split(double x, double& x_hi, double& x_low) {
   double c     = exact_fp_splitter * x;
   double x_big = c - x;
   x_hi  = c - x_big;
   x_low = x - x_hi;
}

/*
 * Fast two-sum.
 * Given |a| > |b|, return nonoverlapping expansion x + y = a + b, where 
 * x is an approximation to a + b and y is the roundoff error.
 */
void exact_fp::fast_two_sum(
   double a, double b, double& x, double& y)
{
   x = a + b;
   double b_v = x - a;
   y = b - b_v;
}

/*
 * Fast two-diff.
 * Given |a| > |b|, return nonoverlapping expansion x + y = a - b, where 
 * x is an approximation to a - b and y is the roundoff error.
 */
void exact_fp::fast_two_diff(
   double a, double b, double& x, double& y)
{
   x = a - b;
   double b_v = a - x;
   y = b_v - b;
}

/*
 * Two-sum.
 * Given a, b, return nonoverlapping expansion x + y = a + b, where
 * x is an approximation to a + b and y is the roundoff error.
 */
void exact_fp::two_sum(
   double a, double b, double& x, double& y)
{
   x = a + b;
   double b_v = x - a;
   double a_v = x - b_v;
   double b_r = b - b_v;
   double a_r = a - a_v;
   y = a_r + b_r;
}

/*
 * Two-diff.
 * Given a, b, return nonoverlapping expansion x + y = a - b, where
 * x is an approximation to a - b and y is the roundoff error.
 */
void exact_fp::two_diff(
   double a, double b, double& x, double& y)
{
   x = a - b;
   double b_v = a - x;
   double a_v = x + b_v;
   double b_r = b_v - b;
   double a_r = a - a_v;
   y = a_r + b_r;
}

/*
 * Two-product.
 * Given a, b, return nonoverlapping expansion x + y = a * b, where
 * x is an approximation to a * b and y is the roundoff error.
 */
void exact_fp::two_product(
   double a, double b, double& x, double& y)
{
   x = a * b;
   double a_hi = 0, a_low = 0;
   double b_hi = 0, b_low = 0;
   exact_fp::split(a, a_hi, a_low);
   exact_fp::split(b, b_hi, b_low);
   double err0 = x - (a_hi * b_hi);
   double err1 = err0 - (a_low * b_hi);
   double err2 = err1 - (a_hi * b_low);
   y = (a_low * b_low) - err2;
}

/*
 * Grow expansion.
 * Given a nonoverlapping expansion and a value, return their sum as a 
 * nonoverlapping expansion.
 */
exact_fp exact_fp::grow_expansion(
   const exact_fp& e, double x)
{
   /* allocate result */
   unsigned long size = e._size + 1;
   exact_fp h(size);
   /* compute sum */
   double Q = x, q = 0, Q_next = 0;
   unsigned long h_pos = 0;
   for (unsigned long n = 0; n < e._size; n++) {
      exact_fp::two_sum(Q, e._fp[n], Q_next, q);
      if (q != 0) { h._fp[h_pos++] = q; }
      Q = Q_next;
   }
   if (Q != 0)
      h._fp[h_pos++] = Q;
   h._size = h_pos;
   return h;
}

/*
 * Scale expansion.
 * Multiply the nonoverlapping expansion by the given value and return the
 * result as a nonoverlapping expansion.
 */
exact_fp exact_fp::scale_expansion(
   const exact_fp& e, double x)
{
   /* check if either argument is zero */
   if ((e._size == 0) || (x == 0))
      return exact_fp();
   /* allocate result */
   unsigned long size = e._size * 2;
   exact_fp h(size);
   /* compute product */
   double Q = 0, q = 0, Q_next = 0;
   double T = 0, t = 0;
   unsigned long h_pos = 0;
   exact_fp::two_product(e._fp[0], x, Q, q);
   if (q != 0)
      h._fp[h_pos++] = q;
   for (unsigned long n = 1; n < e._size; n++) {
      exact_fp::two_product(e._fp[n], x, T, t);
      exact_fp::two_sum(Q, t, Q_next, q);
      if (q != 0) { h._fp[h_pos++] = q; }
      exact_fp::fast_two_sum(T, Q_next, Q, q);
      if (q != 0) { h._fp[h_pos++] = q; }
   }
   if (Q != 0)
      h._fp[h_pos++] = Q;
   h._size = h_pos;
   return h;
}

/*
 * Linear expansion sum.
 * Given two nonoverlapping expansions, return their sum as a nonoverlapping
 * expansion.
 */
exact_fp exact_fp::linear_expansion_sum(
   const exact_fp& e0, const exact_fp& e1)
{
   /* check if either argument is zero */
   if (e0._size == 0)
      return e1;
   if (e1._size == 0)
      return e0;
   /* compute size of resulting expansion */
   unsigned long size = e0._size + e1._size;
   /* merge e0 and e1 into a single ordered sequence g */
   exact_fp g(size);
   unsigned long e0_pos = 0;
   unsigned long e1_pos = 0;
   for (unsigned long n = 0; n < size; n++) {
      if (e0_pos == e0._size) {
         g._fp[n] = e1._fp[e1_pos++];
      } else if (e1_pos == e1._size) {
         g._fp[n] = e0._fp[e0_pos++];
      } else {
         double e0_val = e0._fp[e0_pos];
         double e1_val = e1._fp[e1_pos];
         if (math::abs(e0_val) < math::abs(e1_val)) {
            g._fp[n] = e0_val;
            e0_pos++;
         } else {
            g._fp[n] = e1_val;
            e1_pos++;
         }
      }
   }
   /* compute sum */
   exact_fp h(size);
   double Q = 0, q = 0, Q_next = 0;
   double R = 0, r = 0;
   unsigned long h_pos = 0;
   exact_fp::fast_two_sum(g._fp[1], g._fp[0], Q, q);
   for (unsigned long n = 2; n < size; n++) {
      exact_fp::fast_two_sum(g._fp[n], q, R, r);
      exact_fp::two_sum(Q, R, Q_next, q);
      Q = Q_next;
      if (r != 0)
         h._fp[h_pos++] = r;
   }
   if (q != 0)
      h._fp[h_pos++] = q;
   if (Q != 0)
      h._fp[h_pos++] = Q;
   h._size = h_pos;
   return h;
}

/*
 * Expansion product.
 * Given two nonoverlapping expansions, return their product as a
 * nonoverlapping expansion.
 */
exact_fp exact_fp::expansion_product(
   const exact_fp& e0, const exact_fp& e1)
{
   /* check if either argument is zero */
   if ((e0._size == 0) || (e1._size == 0))
      return exact_fp();
   /* multiply the longer sequence by the shorter one */
   return (
      (e0._size < e1._size) ?
         exact_fp::expansion_product(e1, e0._fp, e0._size)
       : exact_fp::expansion_product(e0, e1._fp, e1._size)
   );
}

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
exact_fp exact_fp::expansion_product(
   const exact_fp& e, double* v, unsigned long n)
{
   /* check if multiplying by scalar */
   if (n == 1) {
      /* scale expansion */
      return exact_fp::scale_expansion(e, v[0]);
   } else {
      /* recursively sum the results of multiplcations */
      unsigned long mid = n / 2;
      exact_fp e0 = exact_fp::expansion_product(e, v,       mid);
      exact_fp e1 = exact_fp::expansion_product(e, &v[mid], n - mid);
      return (e0 + e1);
   }
}

} /* namespace math */
