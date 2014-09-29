/*
 * Math.
 */
#include "lang/types/type_ranges.hh"
#include "math/math.hh"

namespace math {

/*
 * Define inverse hyperbolic trig functions.
 * (if the standard library does not define them)
 */
#if MATH__MATH__DEFINE_ACOSH
float acosh(float x) {
   return log(x + sqrt(x*x - 1));
}

double acosh(double x) {
   return log(x + sqrt(x*x - 1));
}

long double acosh(long double x) {
   return log(x + sqrt(x*x - 1));
}
#endif

#if MATH__MATH__DEFINE_ASINH
float asinh(float x) {
   return log(x + sqrt(x*x + 1));
}

double asinh(double x) {
   return log(x + sqrt(x*x + 1));
}

long double asinh(long double x) {
   return log(x + sqrt(x*x + 1));
}
#endif

#if MATH__MATH__DEFINE_ATANH
float atanh(float x) {
   return (log((1 + x)/(1 - x)))/2;
}

double atanh(double x) {
   return (log((1 + x)/(1 - x)))/2;
}

long double atanh(long double x) {
   return (log((1 + x)/(1 - x)))/2;
}
#endif

#if MATH__MATH__DEFINE_ROUND
float round(float x) {
   return ((x < 0) ? ceil(x - 0.5) : floor(x + 0.5));
}

double round(double x) {
   return ((x < 0) ? ceil(x - 0.5) : floor(x + 0.5));
}

long double round(long double x) {
   return ((x < 0) ? ceil(x - 0.5) : floor(x + 0.5));
}
#endif

/*
 * Return epsilon (numeric precision) for the given built-in numeric type.
 */
template <>
float eps<float>() {
   return FLT_EPSILON;
}

template <>
double eps<double>() {
   return DBL_EPSILON;
}

template <>
long double eps<long double>() {
   return LDBL_EPSILON;
}

} /* namespace math */
