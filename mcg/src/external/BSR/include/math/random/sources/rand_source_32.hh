/*
 * Pseudorandom source for 32-bit generators.
 * This class implements the transformations to all other built-in types for 
 * generators which produce random data in 32-bit (unsigned long) chunks.
 */
#ifndef MATH__RANDOM__SOURCES__RAND_SOURCE_32_HH
#define MATH__RANDOM__SOURCES__RAND_SOURCE_32_HH

#include "math/random/sources/rand_source.hh"

namespace math {
namespace random {
namespace sources {

/*
 * Abstract base class for pseudorandom 32-bit source.
 */
class rand_source_32 : public rand_source {
public:
   /*
    * Destructor.
    */
   virtual ~rand_source_32() = 0;

   /*
    * Generate a pseudorandom boolean uniformly distributed on {true, false}.
    */
   bool gen_uniform_bool();
   
   /*
    * Generate a pseudorandom integer uniformly distributed on the closed
    * interval [0, 2^N - 1], where N is the number of bits of precision in
    * the corresponding datatype.
    */
   unsigned char      gen_uniform_unsigned_char();
   unsigned short     gen_uniform_unsigned_short();
   unsigned int       gen_uniform_unsigned_int();
   unsigned long      gen_uniform_unsigned_long();
   unsigned long long gen_uniform_unsigned_long_long();
   
protected:
   /*
    * Generate a pseudorandom uniformly distributed 32-bit chunk.
    */
   virtual unsigned int generate() = 0;
};

} /* namespace sources */
} /* namespace random */
} /* namespace math */

#endif
