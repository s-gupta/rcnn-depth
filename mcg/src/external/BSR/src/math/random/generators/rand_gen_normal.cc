/*
 * This file implements helper functions for generators that produce 
 * normally distributed pseudorandom numbers.
 */
#include "math/random/generators/rand_gen_normal.hh"
#include "math/random/sources/rand_source.hh"

namespace math {
namespace random {
namespace generators {
/*
 * Imports.
 */
using math::random::sources::rand_source;

/*
 * Pure virtual destructor.
 */
rand_gen_normal_base::~rand_gen_normal_base() { }

/*
 * Overloaded functions for generating normally distributed pseudorandom 
 * elements of the built-in floating point types.
 */
float rand_gen_normal_base::generate(
   rand_source& r,
   float mean, 
   float variance) 
{
   return r.gen_normal_float(mean, variance);
}

double rand_gen_normal_base::generate(
   rand_source& r,
   double mean,
   double variance)
{
   return r.gen_normal_double(mean, variance);
}

long double rand_gen_normal_base::generate(
   rand_source& r,
   long double mean,
   long double variance)
{
   return r.gen_normal_long_double(mean, variance);
}

} /* namespace generators */
} /* namespace random */
} /* namespace math */
