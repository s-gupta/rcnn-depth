/*
 * Functions for random permutations.
 */
#include "collections/abstract/collection.hh"
#include "lang/array.hh"
#include "math/random/generators/rand_gen.hh"
#include "math/random/generators/rand_gen_uniform.hh"
#include "math/random/util/randperm.hh"
#include "math/random/util/sampler.hh"

namespace math {
namespace random {
namespace util {
/*
 * Imports.
 */
using collections::abstract::collection;
using lang::array;
using math::random::generators::rand_gen;
using math::random::generators::rand_gen_uniform;

/*
 * Random permutation.
 * Return a random permutation of the integers {0,1,...,N-1}.
 */
array<unsigned long> randperm(unsigned long N) {
   rand_gen_uniform<> r;
   return randperm(N, r);
}

/*
 * Random permutation.
 * Return a random permutation of the integers {0,1,...,N-1}.
 * Use the given random number generator as the source of randomness.
 */
array<unsigned long> randperm(unsigned long N, rand_gen<>& r) {
   array<double> rand_arr(N);
   for (unsigned long n = 0; n < N; n++)
      rand_arr[n] = r.generate();
   return rand_arr.sort_idx();
}

/*
 * Weighted random permutation.
 *
 * Return a weighted random permutation of the integers {0,1,...,N-1}, where 
 * N is the number of weights in the given collection.  Each integer is put 
 * into a bag and samples are drawn, without replacement, with probability 
 * proportional to the corresponding weight.
 */
array<unsigned long> randperm(const collection<double>& weights) {
   /* create collection of integers to sample */
   unsigned long N = weights.size();
   array<unsigned long> a(N);
   list<unsigned long> items;
   for (unsigned long n = 0; n < N; n++) {
      a[n] = n;
      items.add(a[n]);
   }
   /* create sampler */
   sampler<unsigned long> s(items, weights, false);
   /* sample items by weight */
   array<unsigned long> a_perm(N);
   for (unsigned long n = 0; n < N; n++)
      a_perm[n] = s.sample();
   return a_perm;
}

} /* namespace util */
} /* namespace random */
} /* namespace math */
