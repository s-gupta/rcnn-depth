/*
 * Basic math functions.
 */
#ifndef MATH__MATH_HH
#define MATH__MATH_HH

#include <cmath>

/*
 * Set this definition to (false) if constants are defined in the included
 * system headers.  Set it to (true) if constants need to be defined below.
 */
#define MATH__MATH__DEFINE_CONSTANTS (false)

/*
 * Set these definitions to (false) if the corresponding function on floating
 * point types is defined in the included system headers and (true) if it 
 * needs to be defined below.
 */
#define MATH__MATH__DEFINE_ACOSH (true)
#define MATH__MATH__DEFINE_ASINH (true)
#define MATH__MATH__DEFINE_ATANH (true)
#define MATH__MATH__DEFINE_ROUND (true)

/*
 * Mathematical constants.
 */
#ifdef MATH__MATH__DEFINE_CONSTANTS
#ifndef M_E
#define M_E		2.7182818284590452354	/* e */
#define M_LOG2E	        1.4426950408889634074	/* log_2 e */
#define M_LOG10E	0.43429448190325182765	/* log_10 e */
#define M_LN2		0.69314718055994530942	/* log_e 2 */
#define M_LN10		2.30258509299404568402	/* log_e 10 */
#define M_PI		3.14159265358979323846	/* pi */
#define M_PI_2		1.57079632679489661923	/* pi/2 */
#define M_PI_4		0.78539816339744830962	/* pi/4 */
#define M_1_PI		0.31830988618379067154	/* 1/pi */
#define M_2_PI		0.63661977236758134308	/* 2/pi */
#define M_2_SQRTPI	1.12837916709551257390	/* 2/sqrt(pi) */
#define M_SQRT2	        1.41421356237309504880	/* sqrt(2) */
#define M_SQRT1_2	0.70710678118654752440	/* 1/sqrt(2) */
#endif
#endif

/*
 * Mathematical constants - 128-bit precision.
 */
#ifdef MATH__MATH__DEFINE_CONSTANTS
#define M_El		2.7182818284590452353602874713526625L  /* e */
#define M_LOG2El	1.4426950408889634073599246810018921L  /* log_2 e */
#define M_LOG10El	0.4342944819032518276511289189166051L  /* log_10 e */
#define M_LN2l		0.6931471805599453094172321214581766L  /* log_e 2 */
#define M_LN10l	        2.3025850929940456840179914546843642L  /* log_e 10 */
#define M_PIl		3.1415926535897932384626433832795029L  /* pi */
#define M_PI_2l	        1.5707963267948966192313216916397514L  /* pi/2 */
#define M_PI_4l	        0.7853981633974483096156608458198757L  /* pi/4 */
#define M_1_PIl	        0.3183098861837906715377675267450287L  /* 1/pi */
#define M_2_PIl	        0.6366197723675813430755350534900574L  /* 2/pi */
#define M_2_SQRTPIl     1.1283791670955125738961589031215452L  /* 2/sqrt(pi) */
#define M_SQRT2l	1.4142135623730950488016887242096981L  /* sqrt(2) */
#define M_SQRT1_2l	0.7071067811865475244008443621048490L  /* 1/sqrt(2) */
#endif

#ifdef __APPLE__
#define M_El		2.7182818284590452353602874713526625L  /* e */
#define M_LOG2El	1.4426950408889634073599246810018921L  /* log_2 e */
#define M_LOG10El	0.4342944819032518276511289189166051L  /* log_10 e */
#define M_LN2l		0.6931471805599453094172321214581766L  /* log_e 2 */
#define M_LN10l	        2.3025850929940456840179914546843642L  /* log_e 10 */
#define M_PIl		3.1415926535897932384626433832795029L  /* pi */
#define M_PI_2l	        1.5707963267948966192313216916397514L  /* pi/2 */
#define M_PI_4l	        0.7853981633974483096156608458198757L  /* pi/4 */
#define M_1_PIl	        0.3183098861837906715377675267450287L  /* 1/pi */
#define M_2_PIl	        0.6366197723675813430755350534900574L  /* 2/pi */
#define M_2_SQRTPIl     1.1283791670955125738961589031215452L  /* 2/sqrt(pi) */
#define M_SQRT2l	1.4142135623730950488016887242096981L  /* sqrt(2) */
#define M_SQRT1_2l	0.7071067811865475244008443621048490L  /* 1/sqrt(2) */
#endif

/*
 * Add the standard set of math functions to the math namespace.
 */
namespace math {
/*
 * Export standard math functions.
 */
using std::cos;
using std::sin;
using std::tan;
using std::cosh;
using std::sinh;
using std::tanh;

using std::acos;
using std::asin;
using std::atan;

#ifdef MATH__MATH__DEFINE_ACOSH
#else
using std::acosh;
#endif

#ifdef MATH__MATH__DEFINE_ASINH 
#else
using std::asinh;
#endif

#ifdef MATH__MATH__DEFINE_ATANH
#else
using std::atanh;
#endif

using std::sqrt;
using std::log;
using std::exp;
using std::pow;

using std::atan2;

using std::abs;

using std::floor;
using std::ceil;

#ifdef MATH__MATH__DEFINE_ROUND
#else
using std::round
#endif

/*
 * Declare inverse hyperbolic trig functions.
 * (if the standard library does not define them)
 */
#ifdef MATH__MATH__DEFINE_ACOSH
float acosh(float);
double acosh(double);
long double acosh(long double);
#endif

#ifdef MATH__MATH__DEFINE_ASINH
float asinh(float);
double asinh(double);
long double asinh(long double);
#endif

#ifdef MATH__MATH__DEFINE_ATANH
float atanh(float);
double atanh(double);
long double atanh(long double);
#endif

#ifdef MATH__MATH__DEFINE_ROUND
float round(float);
double round(double);
long double round(long double);
#endif

/*
 * Return epsilon (numeric precision) for the given numeric type.
 *
 * For floating point types, epsilon is the distance from 1.0 to the next
 * largest value of that type.
 *
 * Any user-defined numeric types requiring a nonzero eps() value should
 * implement it through template specialization.
 */
template <typename T> T eps();

template <typename T>
T eps() {
   return T();
}

template <> float eps<float>();
template <> double eps<double>();
template <> long double eps<long double>();

} /* namespace math */

#endif
