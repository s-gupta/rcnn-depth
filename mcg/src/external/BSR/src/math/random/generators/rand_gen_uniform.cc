/*
 * This file implements helper functions for generators that produce 
 * uniformly distributed pseudorandom numbers.
 */
#include "math/random/generators/rand_gen_uniform.hh"
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
rand_gen_uniform_base::~rand_gen_uniform_base() { }

/*
 * Overloaded functions for generating uniformly distributed pseudorandom 
 * elements of the built-in numeric types.
 */
bool rand_gen_uniform_base::generate(rand_source& r, bool min, bool max) {
   return (min == max) ? min : r.gen_uniform_bool();
}

char rand_gen_uniform_base::generate(rand_source& r, char min, char max) {
   return (r.gen_uniform_char() % (max - min + 1)) + min;
}

short rand_gen_uniform_base::generate(rand_source& r, short min, short max) {
   return (r.gen_uniform_short() % (max - min + 1)) + min;
}

int rand_gen_uniform_base::generate(rand_source& r, int min, int max) {
   return (r.gen_uniform_int() % (max - min + 1)) + min;
}

long rand_gen_uniform_base::generate(rand_source& r, long min, long max) {
   return (r.gen_uniform_long() % (max - min + 1)) + min;
}

long long rand_gen_uniform_base::generate(rand_source& r, long long min, long long max) {
   return (r.gen_uniform_long_long() % (max - min + 1)) + min;
}

unsigned char rand_gen_uniform_base::generate(rand_source& r, unsigned char min, unsigned char max) {
   return (r.gen_uniform_char() % (max - min + 1)) + min;
}

unsigned short rand_gen_uniform_base::generate(rand_source& r, unsigned short min, unsigned short max) {
   return (r.gen_uniform_short() % (max - min + 1)) + min;
}

unsigned int rand_gen_uniform_base::generate(rand_source& r, unsigned int min, unsigned int max) {
   return (r.gen_uniform_int() % (max - min + 1)) + min;
}

unsigned long rand_gen_uniform_base::generate(rand_source& r, unsigned long min, unsigned long max) {
   return (r.gen_uniform_long() % (max - min + 1)) + min;
}

unsigned long long rand_gen_uniform_base::generate(rand_source& r, unsigned long long min, unsigned long long max) {
   return (r.gen_uniform_long_long() % (max - min + 1)) + min;
}

float rand_gen_uniform_base::generate(rand_source& r, float min, float max) {
   return (r.gen_uniform_closed_float() * (max - min)) + min;
}

double rand_gen_uniform_base::generate(rand_source& r, double min, double max) {
   return (r.gen_uniform_closed_double() * (max - min)) + min;
}

long double rand_gen_uniform_base::generate(rand_source& r, long double min, long double max) {
   return (r.gen_uniform_closed_long_double() * (max - min)) + min;
}

} /* namespace generators */
} /* namespace random */
} /* namespace math */
