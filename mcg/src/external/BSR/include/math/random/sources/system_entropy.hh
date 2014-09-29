/*
 * Pseudorandom source that uses system entropy to generate random numbers.
 * This generator is useful for initializing other generators with random 
 * seeds.
 */
#ifndef MATH__RANDOM__SOURCES__SYSTEM_ENTROPY_HH
#define MATH__RANDOM__SOURCES__SYSTEM_ENTROPY_HH

#include "math/random/sources/rand_source.hh"

namespace math {
namespace random {
namespace sources {

/*
 * Pseudorandom source that uses system entropy.
 */
class system_entropy : public rand_source {
public:
   /*
    * Constructor.
    */
   system_entropy();

   /*
    * Copy constructor.
    */
   system_entropy(const system_entropy&);
   
   /*
    * Destructor.
    */
   virtual ~system_entropy();

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
};

} /* namespace sources */
} /* namespace random */
} /* namespace math */

#endif
