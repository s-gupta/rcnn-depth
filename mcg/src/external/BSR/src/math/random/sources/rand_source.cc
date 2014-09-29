/*
 * Pseudorandom source.
 */
#include "lang/types/type_ranges.hh"
#include "lang/types/type_sizes.hh"
#include "math/math.hh"
#include "math/random/sources/rand_source.hh"

namespace math {
namespace random {
namespace sources {

/*
 * Pure virtual destructor.
 */
rand_source::~rand_source() { }

/*
 * Generate a pseudorandom integer uniformly distributed on the closed
 * interval [0, 2^(N-1) - 1], where N is the number of bits of precision in
 * the corresponding unsigned datatype.
 */
char rand_source::gen_uniform_char() {
   return static_cast<char>(this->gen_uniform_unsigned_char() >> 1);
}

short rand_source::gen_uniform_short() {
   return static_cast<short>(this->gen_uniform_unsigned_short() >> 1);
}

int rand_source::gen_uniform_int() {
   return static_cast<int>(this->gen_uniform_unsigned_int() >> 1);
}

long rand_source::gen_uniform_long() {
   return static_cast<long>(this->gen_uniform_unsigned_long() >> 1);
}

long long rand_source::gen_uniform_long_long() {
   return static_cast<long long>(this->gen_uniform_unsigned_long_long() >> 1);
}

/*
 * Generate a pseudorandom real number uniformly distributed on the 
 * closed interval [0,1].
 */
float rand_source::gen_uniform_closed_float() {
   return (this->gen_uniform_unsigned_int() >> (UINT_BITS - FLT_MANT_DIG)) * 
          1.0/(static_cast<float>(UINT_MAX >> (UINT_BITS - FLT_MANT_DIG)));
}

double rand_source::gen_uniform_closed_double() {
   return (this->gen_uniform_unsigned_long_long() >> (ULONG_LONG_BITS - DBL_MANT_DIG)) * 
          1.0/(static_cast<double>(ULONG_LONG_MAX >> (ULONG_LONG_BITS - DBL_MANT_DIG)));
}

long double rand_source::gen_uniform_closed_long_double() {
   return static_cast<long double>(this->gen_uniform_closed_double());
}

/*
 * Generate a pseudorandom real number uniformly distributed on the
 * half-open interval [0,1).
 */
float rand_source::gen_uniform_half_open_float() {
   return (this->gen_uniform_unsigned_int() >> (UINT_BITS - FLT_MANT_DIG)) * 
          1.0/(static_cast<float>(UINT_MAX >> (UINT_BITS - FLT_MANT_DIG)) + 1.0);
}

double rand_source::gen_uniform_half_open_double() {
   return (this->gen_uniform_unsigned_long_long() >> (ULONG_LONG_BITS - DBL_MANT_DIG)) * 
          1.0/(static_cast<double>(ULONG_LONG_MAX >> (ULONG_LONG_BITS - DBL_MANT_DIG)) + 1.0);
}

long double rand_source::gen_uniform_half_open_long_double() {
   return static_cast<long double>(this->gen_uniform_half_open_double());
}

/*
 * Generate a pseudorandom real number uniformly distributed on the
 * open interval (0,1).
 */
float rand_source::gen_uniform_open_float() {
   return (static_cast<float>(this->gen_uniform_unsigned_int() >> (UINT_BITS - FLT_MANT_DIG + 1)) + 0.5) * 
          1.0/(static_cast<float>(UINT_MAX >> (UINT_BITS - FLT_MANT_DIG + 1)) + 1.0);
}

double rand_source::gen_uniform_open_double() {
   return (static_cast<double>(this->gen_uniform_unsigned_long_long() >> (ULONG_LONG_BITS - DBL_MANT_DIG + 1)) + 0.5) * 
          1.0/(static_cast<double>(ULONG_LONG_MAX >> (ULONG_LONG_BITS - DBL_MANT_DIG + 1)) + 1.0);
}

long double rand_source::gen_uniform_open_long_double() {
   return static_cast<long double>(this->gen_uniform_open_double());
}

/*
 * Generate a pseudorandom float normally distributed with mean zero and 
 * unit variance.  Use the Box-Muller transformation.
 */
float rand_source::gen_normal_float() {
   float x, y, r;
   do {
      /* generate x and y uniformly distributed in [-1,1] */
      x = (float(2.0))*(this->gen_uniform_closed_float()) - (float(1.0));
      y = (float(2.0))*(this->gen_uniform_closed_float()) - (float(1.0));
      /* compute r */
      r = (x*x) + (y*y);
   } while ((r == float(0.0)) || (r > float(1.0)));
   return (x * math::sqrt(float(-2.0)*(math::log(r))/r));
}

/*
 * Generate a pseudorandom double normally distributed with mean zero and 
 * unit variance.  Use the Box-Muller transformation.
 */
double rand_source::gen_normal_double() {
   double x, y, r;
   do {
      /* generate x and y uniformly distributed in [-1,1] */
      x = (double(2.0))*(this->gen_uniform_closed_double()) - (double(1.0));
      y = (double(2.0))*(this->gen_uniform_closed_double()) - (double(1.0));
      /* compute r */
      r = (x*x) + (y*y);
   } while ((r == double(0.0)) || (r > double(1.0)));
   return (x * math::sqrt(double(-2.0)*(math::log(r))/r));
}

/*
 * Generate a pseudorandom long double normally distributed with mean zero and 
 * unit variance.  Use the Box-Muller transformation.
 */
long double rand_source::gen_normal_long_double() {
   long double x, y, r;
   do {
      /* generate x and y uniformly distributed in [-1,1] */
      x = (static_cast<long double>(2.0))*(this->gen_uniform_closed_long_double()) - (static_cast<long double>(1.0));
      y = (static_cast<long double>(2.0))*(this->gen_uniform_closed_long_double()) - (static_cast<long double>(1.0));
      /* compute r */
      r = (x*x) + (y*y);
   } while ((r == static_cast<long double>(0.0)) || (r > static_cast<long double>(1.0)));
   return (x * math::sqrt((static_cast<long double>(-2.0))*(math::log(r))/r));
}

/*
 * Generate a pseudorandom float normally distributed 
 * with the specified mean and variance.
 */
float rand_source::gen_normal_float(float mean, float variance) {
   return (mean + (math::sqrt(variance))*(this->gen_normal_float()));
}

/*
 * Generate a pseudorandom double normally distributed 
 * with the specified mean and variance.
 */
double rand_source::gen_normal_double(double mean, double variance) {
   return (mean + (math::sqrt(variance))*(this->gen_normal_double()));
}

/*
 * Generate a pseudorandom long double normally distributed 
 * with the specified mean and variance.
 */
long double rand_source::gen_normal_long_double(long double mean, long double variance) {
   return (mean + (math::sqrt(variance))*(this->gen_normal_long_double()));
}

} /* namespace sources */
} /* namespace random */
} /* namespace math */
