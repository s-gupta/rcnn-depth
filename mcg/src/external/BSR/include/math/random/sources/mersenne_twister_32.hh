/*
 * Mersenne Twister 32-bit pseudorandom source.
 * 
 * The Mersenne Twister is a fast pseudorandom number generator with a period 
 * of 2^19937-1.  It is not cryptographically secure as future outputs can be 
 * predicted given a sufficiently long subsequence of past outputs.
 */
#ifndef MATH__RANDOM__SOURCES__MERSENNE_TWISTER_32_HH
#define MATH__RANDOM__SOURCES__MERSENNE_TWISTER_32_HH

#include "lang/array.hh"
#include "math/random/sources/rand_source.hh"
#include "math/random/sources/rand_source_32.hh"

namespace math {
namespace random {
namespace sources {
/*
 * Imports.
 */
using lang::array;

/*
 * Mersenne Twister pseudorandom source.
 */
class mersenne_twister_32 : public rand_source_32 {
public:
   /*
    * Constructor.
    * Initialize the state randomly using entropy from the system.
    */
   mersenne_twister_32();

   /*
    * Constructor.
    * Initialize the state by using the given random source.
    */
   explicit mersenne_twister_32(rand_source&);
   
   /*
    * Constructor.
    * Initialize the state using the given data.
    */
   explicit mersenne_twister_32(const array<unsigned int>&);
   
   /*
    * Copy constructor.
    * Give the copy the same state as the original.
    */
   mersenne_twister_32(const mersenne_twister_32&);

   /*
    * Destructor.
    */
   virtual ~mersenne_twister_32();

protected:
   /*
    * Generator state.
    */
   array<unsigned int> _state;   /* random state */
   unsigned long _curr;          /* current position within state */

   /*
    * Initialize the state using the given random source or given data.
    */
   void initialize(rand_source&);
   void initialize(const array<unsigned int>&);

   /*
    * Generate a pseudorandom unsigned 32-bit integer uniformly distributed 
    * on the closed interval [0, 2^32-1].
    */
   unsigned int generate();
};

} /* namespace sources */
} /* namespace random */
} /* namespace math */

#endif
