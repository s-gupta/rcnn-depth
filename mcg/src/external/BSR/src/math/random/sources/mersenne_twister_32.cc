/*
 * Mersenne Twister 32-bit pseudorandom source.
 * 
 * This implementation is based on freely distributed code by Takuji Nishimura
 * and Makoto Matsumoto.  The copyright notice for the original version is 
 * included below.
 */

/* 
 * A C-program for MT19937, with initialization improved 2002/1/26.
 * Coded by Takuji Nishimura and Makoto Matsumoto.
 *
 * Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
 * All rights reserved.                          
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *   1. Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *
 *   3. The names of its contributors may not be used to endorse or promote 
 *      products derived from this software without specific prior written 
 *      permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Any feedback is very welcome.
 * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
 * email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
 */

#include "lang/array.hh"
#include "math/random/sources/mersenne_twister_32.hh"
#include "math/random/sources/rand_source.hh"
#include "math/random/sources/system_entropy.hh"

namespace math {
namespace random {
namespace sources {
/*
 * Imports.
 */
using lang::array;

/*
 * Generator parameters.
 */
namespace {
const unsigned long N = 624;
const unsigned long M = 397;
const unsigned int MATRIX_A = 0x9908b0df;    /* constant vector a */
const unsigned int UPPER_MASK = 0x80000000;  /* most significant w-r bits */
const unsigned int LOWER_MASK = 0x7fffffff;  /* least significant r bits */
}

/*
 * Constructor.
 * Initialize the state randomly using entropy from the system.
 */
mersenne_twister_32::mersenne_twister_32()
 : _state(N),
   _curr(N)
{
   system_entropy sys_entropy;
   this->initialize(sys_entropy);
}

/*
 * Constructor.
 * Initialize the state by using the given random source.
 */
mersenne_twister_32::mersenne_twister_32(rand_source& r)
 : _state(N),
   _curr(N)
{
   this->initialize(r);
}

/*
 * Constructor.
 * Initialize the state using the given data.
 */
mersenne_twister_32::mersenne_twister_32(const array<unsigned int>& seed)
 : _state(N),
   _curr(N)
{
   this->initialize(seed);
}

/*
 * Initialize the state using the given random source.
 */
void mersenne_twister_32::initialize(rand_source& r) {
   /* generate random data for initialization */
   array<unsigned int> seed(N);
   for (unsigned long n = 0; n < N; n++)
      seed[n] = r.gen_uniform_unsigned_int();
   /* initialize using random data */
   this->initialize(seed);
}

/*
 * Initialize the state using the given data.
 */
void mersenne_twister_32::initialize(const array<unsigned int>& seed) {
   /* check if seed is empty */
   if (seed.is_empty()) {
      /* use a predefined seed */
      array<unsigned int> default_seed(4);
      default_seed[0] = 0x123;
      default_seed[1] = 0x234;
      default_seed[2] = 0x345;
      default_seed[3] = 0x456;
      this->initialize(default_seed);
   } else {
      /* initialize state to predefined pattern */
      _state[0]= 19650218;
      for (unsigned long n = 1; n < N; n++) {
         _state[n] = 
            (1812433253 * (_state[n-1] ^ (_state[n-1] >> 30)) + n); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array _state[].                    */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        _state[n] &= 0xffffffff;
        /* for >32 bit machines */
      }
      /* modify state with the given seed */
      unsigned long seed_size = seed.size();
      unsigned long i = 1;
      unsigned long j = 0;
      unsigned long k = (N > seed_size) ? N : seed_size;
      for (; k > 0; k--) {
         _state[i] = (_state[i] ^ ((_state[i-1] ^ (_state[i-1] >> 30)) * 1664525))
            + seed[j] + static_cast<unsigned int>(j); /* non linear */
         _state[i] &= 0xffffffff; /* for WORDSIZE > 32 machines */
         i++; j++;
         if (i >= N) { _state[0] = _state[N-1]; i = 1; }
         if (j >= seed_size) { j = 0; }
      }
      for (k = N-1; k > 0; k--) {
         _state[i] = (_state[i] ^ ((_state[i-1] ^ (_state[i-1] >> 30)) * 1566083941))
            - static_cast<unsigned int>(i); /* non linear */
         _state[i] &= 0xffffffff; /* for WORDSIZE > 32 machines */
         i++;
         if (i >= N) { _state[0] = _state[N-1]; i = 1; }
      }
      _state[0] = 0x80000000; /* MSB is 1; assuring non-zero initial array */ 
   }
}

/*
 * Copy constructor.
 * Give the copy the same state as the original.
 */
mersenne_twister_32::mersenne_twister_32(const mersenne_twister_32& r)
 : _state(r._state),
   _curr(r._curr)
{ }

/*
 * Destructor.
 */
mersenne_twister_32::~mersenne_twister_32() {
   /* do nothing */
}

/*
 * Generate a pseudorandom unsigned 32-bit integer uniformly distributed 
 * on the closed interval [0, 2^32-1].
 */
unsigned int mersenne_twister_32::generate() {
   static const unsigned int mag01[2] = {0x0, MATRIX_A};
   unsigned int x;
   /* check if new items need to be generated */
   if (_curr >= N) {
      /* generate NN items at once */
      unsigned long i;
      for (i = 0; i < (N-M); i++) {
         x = (_state[i] & UPPER_MASK) | (_state[i+1] & LOWER_MASK);
         _state[i] = _state[i+M] ^ (x >> 1) ^ mag01[static_cast<unsigned long>(x & 0x1)];
      }
      for (; i < (N-1); i++) {
         x = (_state[i] & UPPER_MASK) | (_state[i+1] & LOWER_MASK);
         _state[i] = _state[i + M - N] ^ (x >> 1) ^ mag01[static_cast<unsigned long>(x & 0x1)];
      }
      x = (_state[N-1] & UPPER_MASK) | (_state[0] & LOWER_MASK);
      _state[N-1] = _state[M-1] ^ (x >> 1) ^ mag01[static_cast<unsigned long>(x & 0x1)];
      _curr = 0;
   }
   /* retrieve next item */
   x = _state[_curr++];
   x ^= (x >> 11);
   x ^= (x <<  7) & 0x9d2c5680;
   x ^= (x << 15) & 0xefc60000;
   x ^= (x >> 18);
   return x;
}

} /* namespace sources */
} /* namespace random */
} /* namespace math */
