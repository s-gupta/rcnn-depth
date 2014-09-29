/*
 * Pseudorandom source for 64-bit generators.
 */
#include "lang/types/type_ranges.hh"
#include "lang/types/type_sizes.hh"
#include "math/random/sources/rand_source_64.hh"

namespace math {
namespace random {
namespace sources {

/*
 * Pure virtual destructor.
 */
rand_source_64::~rand_source_64() { }

/*
 * Generate a pseudorandom boolean uniformly distributed on {true, false}.
 */
bool rand_source_64::gen_uniform_bool() {
   return static_cast<bool>(this->generate() >> (ULONG_LONG_BITS - 1));
}

/*
 * Generate a pseudorandom integer uniformly distributed on the closed
 * interval [0, 2^N - 1], where N is the number of bits of precision in
 * the corresponding datatype.
 */
unsigned char rand_source_64::gen_uniform_unsigned_char() {
   return static_cast<unsigned char>(this->generate() >> (ULONG_LONG_BITS - UCHAR_BITS));
}

unsigned short rand_source_64::gen_uniform_unsigned_short() {
   return static_cast<unsigned short>(this->generate() >> (ULONG_LONG_BITS - USHRT_BITS));
}

unsigned int rand_source_64::gen_uniform_unsigned_int() {
   return static_cast<unsigned int>(this->generate() >> (ULONG_LONG_BITS - UINT_BITS));
}

unsigned long rand_source_64::gen_uniform_unsigned_long() {
   return static_cast<unsigned long>(this->generate() >> (ULONG_LONG_BITS - ULONG_BITS));
}

unsigned long long rand_source_64::gen_uniform_unsigned_long_long() {
   return this->generate();
}

} /* namespace sources */
} /* namespace random */
} /* namespace math */
