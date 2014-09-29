/*
 * Pseudorandom source for 64-bit generators.
 * This class implements the transformations to all other built-in types for 
 * generators which produce random data in 64-bit (unsigned long long) chunks.
 */
#ifndef MATH__RANDOM__SOURCES__RAND_SOURCE_64_HH
#define MATH__RANDOM__SOURCES__RAND_SOURCE_64_HH

#include "math/random/sources/rand_source.hh"

namespace math {
namespace random {
namespace sources {

/*
 * Abstract base class for pseudorandom 64-bit source.
 */
class rand_source_64 : public rand_source {
public:
   /*
    * Destructor.
    */
   virtual ~rand_source_64() = 0;

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
    * Generate a pseudorandom uniformly distributed 64-bit chunk.
    */
   virtual unsigned long long generate() = 0;
};

} /* namespace sources */
} /* namespace random */
} /* namespace math */

#endif
