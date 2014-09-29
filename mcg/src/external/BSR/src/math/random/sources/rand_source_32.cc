/*
 * Pseudorandom source for 32-bit generators.
 */
#include "lang/types/type_ranges.hh"
#include "lang/types/type_sizes.hh"
#include "math/random/sources/rand_source_32.hh"

namespace math {
namespace random {
namespace sources {

/*
 * Pure virtual destructor.
 */
rand_source_32::~rand_source_32() { }

/*
 * Generate a pseudorandom boolean uniformly distributed on {true, false}.
 */
bool rand_source_32::gen_uniform_bool() {
   return static_cast<bool>(this->generate() >> (UINT_BITS - 1));
}

/*
 * Generate a pseudorandom integer uniformly distributed on the closed
 * interval [0, 2^N - 1], where N is the number of bits of precision in
 * the corresponding datatype.
 */
unsigned char rand_source_32::gen_uniform_unsigned_char() {
   return static_cast<unsigned char>(this->generate() >> (UINT_BITS - UCHAR_BITS));
}

unsigned short rand_source_32::gen_uniform_unsigned_short() {
   return static_cast<unsigned short>(this->generate() >> (UINT_BITS - USHRT_BITS));
}

unsigned int rand_source_32::gen_uniform_unsigned_int() {
   return this->generate();
}

unsigned long rand_source_32::gen_uniform_unsigned_long() {
   /* unsigned long could be either 32 or 64 bits */
   unsigned long long r = this->gen_uniform_unsigned_long_long();
   return static_cast<unsigned long>(r >> (ULONG_LONG_BITS - UINT_BITS));
}

unsigned long long rand_source_32::gen_uniform_unsigned_long_long() {
   /* generate most significant and least significant chunks */
   unsigned long long ms_chunk = static_cast<unsigned long long>(this->generate());
   unsigned long long ls_chunk = static_cast<unsigned long long>(this->generate());
   /* combine them */
   return (ms_chunk << UINT_BITS) + ls_chunk;
}

} /* namespace sources */
} /* namespace random */
} /* namespace math */
