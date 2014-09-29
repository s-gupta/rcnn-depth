/*
 * Default pseudorandom source.
 * This is the default pseudorandom number source used by the generators 
 * found in the math::random::generators namespace.
 */
#ifndef MATH__RANDOM__SOURCES__RAND_SOURCE_DEFAULT_HH
#define MATH__RANDOM__SOURCES__RAND_SOURCE_DEFAULT_HH

#include "lang/pointers/auto_ptr.hh"
#include "math/random/sources/rand_source.hh"

namespace math {
namespace random {
namespace sources {
/*
 * Impotrs.
 */
using lang::pointers::auto_ptr;

/*
 * Default pseudorandom source.
 */
class rand_source_default : public rand_source {
public:
   /*
    * Constructor.
    * Create a randomly initialized pseudorandom source of the default type.
    */
   rand_source_default();

   /*
    * Copy constructor.
    */
   rand_source_default(const rand_source_default&);

   /*
    * Destructor.
    */
   virtual ~rand_source_default();

   /*
    * Generate a pseudorandom boolean uniformly distributed on {true, false}.
    */
   bool gen_uniform_bool();
   
   /*
    * Generate a pseudorandom integer uniformly distributed on the closed
    * interval [0, 2^(N-1) - 1], where N is the number of bits of precision in
    * the corresponding unsigned datatype.
    */
   char      gen_uniform_char();
   short     gen_uniform_short();
   int       gen_uniform_int();
   long      gen_uniform_long();
   long long gen_uniform_long_long();

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

   /*
    * Generate a pseudorandom real number uniformly distributed on the 
    * closed interval [0,1].
    */
   float       gen_uniform_closed_float();
   double      gen_uniform_closed_double();
   long double gen_uniform_closed_long_double();
   
   /*
    * Generate a pseudorandom real number uniformly distributed on the
    * half-open interval [0,1).
    */
   float       gen_uniform_half_open_float();
   double      gen_uniform_half_open_double();
   long double gen_uniform_half_open_long_double();

   /*
    * Generate a pseudorandom real number uniformly distributed on the
    * open interval (0,1).
    */
   float       gen_uniform_open_float();
   double      gen_uniform_open_double();
   long double gen_uniform_open_long_double();

protected:
   auto_ptr<rand_source> _r;  /* random source to use */
};

} /* namespace sources */
} /* namespace random */
} /* namespace math */

#endif
