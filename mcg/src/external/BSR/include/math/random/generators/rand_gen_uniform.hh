/*
 * Generator that produces uniformly distributed pseudorandom numbers.
 */
#ifndef MATH__RANDOM__GENERATORS__RAND_GEN_UNIFORM_HH
#define MATH__RANDOM__GENERATORS__RAND_GEN_UNIFORM_HH

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
 * Base class containing code common to all rand_gen_uniform<T,Syn> templates.
 */
class rand_gen_uniform_base {
public:
   /*
    * Destructor.
    */
   virtual ~rand_gen_uniform_base() = 0;

protected:
   /*
    * Overloaded functions for generating uniformly distributed pseudorandom 
    * elements of the built-in numeric types.
    */
   static bool               generate(rand_source&, bool,               bool);

   static char               generate(rand_source&, char,               char);
   static short              generate(rand_source&, short,              short);
   static int                generate(rand_source&, int,                int);
   static long               generate(rand_source&, long,               long);
   static long long          generate(rand_source&, long long,          long long);

   static unsigned char      generate(rand_source&, unsigned char,      unsigned char);
   static unsigned short     generate(rand_source&, unsigned short,     unsigned short);
   static unsigned int       generate(rand_source&, unsigned int,       unsigned int);
   static unsigned long      generate(rand_source&, unsigned long,      unsigned long);
   static unsigned long long generate(rand_source&, unsigned long long, unsigned long long);

   static float              generate(rand_source&, float,              float);
   static double             generate(rand_source&, double,             double);
   static long double        generate(rand_source&, long double,        long double);
};

/*
 * Generator that produces uniformly distributed pseudorandom numbers.
 *
 * This generator works for the built-in numeric types.
 * Template specialization should be used to implement generators for 
 * other types.
 */
template <typename T = double, typename Syn = unsynchronized>
class rand_gen_uniform : public rand_gen<T>,
                         protected rand_gen_uniform_base,
                         protected Syn {
public:
   /*
    * Create a generator that produces uniformly distributed numbers in the 
    * range given by [min, max] (inclusive).  The specified pseudorandom 
    * source is used as the source of randomness.
    */
   explicit rand_gen_uniform(
      T = T(0),   /* min */
      T = T(1),   /* max */
      auto_ptr<rand_source> = auto_ptr<rand_source>(new rand_source_default())
   );

   /*
    * Copy constructor.
    * Create a generator with the same range as the original.
    */
   rand_gen_uniform(
      const rand_gen_uniform<T,Syn>&,
      auto_ptr<rand_source> = auto_ptr<rand_source>(new rand_source_default())
   );

   /*
    * Destructor.
    */
   virtual ~rand_gen_uniform();

   /*
    * Generate a pseudorandom number.
    */
   T generate();

protected:
   /*
    * Generator parameters.
    */
   T _min;                    /* lower bound of range */
   T _max;                    /* upper bound of range */
   auto_ptr<rand_source> _r;  /* source of randomness */
};

/***************************************************************************
 * Implementation of rand_gen_uniform.
 ***************************************************************************/

/*
 * Create a generator that produces uniformly distributed numbers in the 
 * range given by [min, max] (inclusive).  The specified pseudorandom 
 * source is used as the source of randomness.
 */
template <typename T, typename Syn>
rand_gen_uniform<T,Syn>::rand_gen_uniform(
   T min,
   T max,
   auto_ptr<rand_source> r)
 : Syn(),
   _min(min),
   _max(max),
   _r(r)
{ }

/*
 * Copy constructor.
 * Create a generator with the same range as the original.
 */
template <typename T, typename Syn>
rand_gen_uniform<T,Syn>::rand_gen_uniform(
   const rand_gen_uniform<T,Syn>& gen,
   auto_ptr<rand_source> r)
 : Syn(),
   _min(gen._min),
   _max(gen._max),
   _r(r)
{ }

/*
 * Destructor.
 */
template <typename T, typename Syn>
rand_gen_uniform<T,Syn>::~rand_gen_uniform() {
   /* do nothing */
}

/*
 * Generate a pseudorandom number.
 */
template <typename T, typename Syn>
T rand_gen_uniform<T,Syn>::generate() {
   auto_lock<const Syn> xlock(*this);
   return rand_gen_uniform_base::generate(*_r, _min, _max);
}

} /* namespace generators */
} /* namespace random */
} /* namespace math */

#endif
