/*
 * Functions for random permutations.
 */
#ifndef MATH__RANDOM__UTIL__RANDPERM_HH
#define MATH__RANDOM__UTIL__RANDPERM_HH

#include "collections/abstract/collection.hh"
#include "lang/array.hh"
#include "math/random/generators/rand_gen.hh"

namespace math {
namespace random {
namespace util {
/*
 * Imports.
 */
using collections::abstract::collection;
using lang::array;
using math::random::generators::rand_gen;

/*
 * Random permutation.
 * Return a random permutation of the integers {0,1,...,N-1}.
 */
array<unsigned long> randperm(unsigned long /* N */);

/*
 * Random permutation.
 * Return a random permutation of the integers {0,1,...,N-1}.
 * Use the given random number generator as the source of randomness.
 */
array<unsigned long> randperm(unsigned long /* N */, rand_gen<>&);

/*
 * Weighted random permutation.
 *
 * Return a weighted random permutation of the integers {0,1,...,N-1}, where 
 * N is the number of weights in the given collection.  Each integer is put 
 * into a bag and samples are drawn, without replacement, with probability 
 * proportional to the corresponding weight.
 */
array<unsigned long> randperm(const collection<double>&);

} /* namespace util */
} /* namespace random */
} /* namespace math */

#endif
