/*
 * Abstract interface for pseudorandom number generators.
 */
#ifndef MATH__RANDOM__GENERATORS__RAND_GEN_HH
#define MATH__RANDOM__GENERATORS__RAND_GEN_HH

namespace math {
namespace random {
namespace generators {

/*
 * Abstract interface for pseudorandom number generator.
 */
template <typename T = double>
class rand_gen {
public:
   /*
    * Destructor.
    */
   virtual ~rand_gen() = 0;

   /*
    * Generate a pseudorandom number.
    */
   virtual T generate() = 0;
};

/*
 * Pure virtual destructor.
 */
template <typename T>
rand_gen<T>::~rand_gen() { }

} /* namespace generators */
} /* namespace random */
} /* namespace math */

#endif
