/*
 * Generator that produces normally distributed pseudorandom numbers.
 */
#ifndef MATH__RANDOM__GENERATORS__RAND_GEN_NORMAL_HH
#define MATH__RANDOM__GENERATORS__RAND_GEN_NORMAL_HH

#include "concurrent/threads/synchronization/locks/auto_lock.hh"
#include "concurrent/threads/synchronization/synchronizables/unsynchronized.hh"
#include "lang/pointers/auto_ptr.hh"
#include "math/random/generators/rand_gen.hh"
#include "math/random/sources/rand_source.hh"
#include "math/random/sources/rand_source_default.hh"

namespace math {
namespace random {
namespace generators {
/*
 * Imports.
 */
using concurrent::threads::synchronization::locks::auto_lock;
using concurrent::threads::synchronization::synchronizables::unsynchronized;
using lang::pointers::auto_ptr;
using math::random::sources::rand_source;
using math::random::sources::rand_source_default;

/*
 * Base class containing code common to all rand_gen_normal<T,Syn> templates.
 */
class rand_gen_normal_base {
public:
   /*
    * Destructor.
    */
   virtual ~rand_gen_normal_base() = 0;

protected:
   /*
    * Overloaded functions for generating normally distributed pseudorandom 
    * elements of the built-in floating point types.
    */
   static float       generate(rand_source&, float,       float);
   static double      generate(rand_source&, double,      double);
   static long double generate(rand_source&, long double, long double);
};

/*
 * Generator that produces normally distributed pseudorandom numbers.
 *
 * This generator works for the built-in floating point types.
 * Template specialization should be used to implement generators for 
 * other types.
 */
template <typename T = double, typename Syn = unsynchronized>
class rand_gen_normal : public rand_gen<T>,
                        protected rand_gen_normal_base,
                        protected Syn {
public:
   /*
    * Create a generator with the given mean and variance.
    * The specified pseudorandom source is used as the source of randomness.
    */
   explicit rand_gen_normal(
      T = T(0),   /* mean */
      T = T(1),   /* variance */
      auto_ptr<rand_source> = auto_ptr<rand_source>(new rand_source_default())
   );
   
   /*
    * Copy constructor.
    * Create a generator with the same mean and variance as the original.
    */
   rand_gen_normal(
      const rand_gen_normal<T,Syn>&,
      auto_ptr<rand_source> = auto_ptr<rand_source>(new rand_source_default())
   );

   /*
    * Destructor.
    */
   virtual ~rand_gen_normal();

   /*
    * Generate a pseudorandom number.
    */
   T generate();

protected:
   /*
    * Generator parameters.
    */
   T _mean;                   /* mean of normal distribution */
   T _variance;               /* variance of normal distribution */
   auto_ptr<rand_source> _r;  /* source of randomness */
};

/***************************************************************************
 * Implementation of rand_gen_normal.
 ***************************************************************************/

/*
 * Create a generator with the given mean and variance.
 * The specified pseudorandom source is used as the source of randomness.
 */
template <typename T, typename Syn>
rand_gen_normal<T,Syn>::rand_gen_normal(
   T mean,
   T variance,
   auto_ptr<rand_source> r)
 : Syn(),
   _mean(mean),
   _variance(variance),
   _r(r)
{ }

/*
 * Copy constructor.
 * Create a generator with the same mean and variance as the original.
 */
template <typename T, typename Syn>
rand_gen_normal<T,Syn>::rand_gen_normal(
   const rand_gen_normal<T,Syn>& gen,
   auto_ptr<rand_source> r)
 : Syn(),
   _mean(gen._mean),
   _variance(gen._variance),
   _r(r)
{ }

/*
 * Destructor.
 */
template <typename T, typename Syn>
rand_gen_normal<T,Syn>::~rand_gen_normal() {
   /* do nothing */
}

/*
 * Generate a pseudorandom number.
 */
template <typename T, typename Syn>
T rand_gen_normal<T,Syn>::generate() {
   auto_lock<const Syn> xlock(*this);
   return rand_gen_normal_base::generate(*_r, _mean, _variance);
}

} /* namespace generators */
} /* namespace random */
} /* namespace math */

#endif
