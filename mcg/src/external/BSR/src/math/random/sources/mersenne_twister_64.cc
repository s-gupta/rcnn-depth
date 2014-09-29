/*
 * Mersenne Twister 64-bit pseudorandom source.
 * 
 * This implementation is based on freely distributed code by Takuji Nishimura
 * and Makoto Matsumoto.  The copyright notice for the original version is 
 * included below.
 */
 
/* 
 * A C-program for MT19937-64 (2004/9/29 version).
 * Coded by Takuji Nishimura and Makoto Matsumoto.
 *
 * This is a 64-bit version of Mersenne Twister pseudorandom number
 * generator.
 *
 * Copyright (C) 2004, Makoto Matsumoto and Takuji Nishimura,
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
 * References:
 * T. Nishimura, ``Tables of 64-bit Mersenne Twisters''
 *   ACM Transactions on Modeling and 
 *   Computer Simulation 10. (2000) 348--357.
 * M. Matsumoto and T. Nishimura,
 *   ``Mersenne Twister: a 623-dimensionally equidistributed
 *     uniform pseudorandom number generator''
 *   ACM Transactions on Modeling and 
 *   Computer Simulation 8. (Jan. 1998) 3--30.
 *
 * Any feedback is very welcome.
 * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
 * email: m-mat @ math.sci.hiroshima-u.ac.jp (remove spaces)
 */
#include "lang/array.hh"
#include "math/random/sources/mersenne_twister_64.hh"
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
const unsigned long NN = 312;   /* length of state vector */
const unsigned long MM = 156;   
const unsigned long long MATRIX_A = 0xB5026F5AA96619E9ULL;
const unsigned long long UM = 0xFFFFFFFF80000000ULL; /* most significant 33 bits */
const unsigned long long LM = 0x7FFFFFFFULL;         /* least significant 31 bits */
}

/*
 * Constructor.
 * Initialize the state randomly using entropy from the system.
 */
mersenne_twister_64::mersenne_twister_64()
 : _state(NN),
   _curr(NN)
{
   system_entropy sys_entropy;
   this->initialize(sys_entropy);
}

/*
 * Constructor.
 * Initialize the state by using the given random source.
 */
mersenne_twister_64::mersenne_twister_64(rand_source& r)
 : _state(NN),
   _curr(NN)
{
   this->initialize(r);
}

/*
 * Constructor.
 * Initialize the state using the given data.
 */
mersenne_twister_64::mersenne_twister_64(const array<unsigned long long>& seed)
 : _state(NN),
   _curr(NN)
{
   this->initialize(seed);
}

/*
 * Initialize the state using the given random source.
 */
void mersenne_twister_64::initialize(rand_source& r) {
   /* generate random data for initialization */
   array<unsigned long long> seed(NN);
   for (unsigned long n = 0; n < NN; n++)
      seed[n] = r.gen_uniform_unsigned_long_long();
   /* initialize using random data */
   this->initialize(seed);
}

/*
 * Initialize the state using the given data.
 */
void mersenne_twister_64::initialize(const array<unsigned long long>& seed) {
   /* check if seed is empty */
   if (seed.is_empty()) {
      /* use a predefined seed */
      array<unsigned long long> default_seed(4);
      default_seed[0] = 0x12345ULL;
      default_seed[1] = 0x23456ULL;
      default_seed[2] = 0x34567ULL;
      default_seed[3] = 0x45678ULL;
      this->initialize(default_seed);
   } else {
      /* initialize state to predefined pattern */
      _state[0] = 19650218ULL;
      for (unsigned long n = 1; n < NN; n++)
        _state[n] = (6364136223846793005ULL * (_state[n-1] ^ (_state[n-1] >> 62)) + n);
      /* modify state with the given seed */
      unsigned long seed_size = seed.size();
      unsigned long i = 1;
      unsigned long j = 0;
      unsigned long k = (NN > seed_size) ? NN : seed_size;
      for (; k > 0; k--) {
         _state[i] = (_state[i] ^ ((_state[i-1] ^ (_state[i-1] >> 62)) * 3935559000370003845ULL))
            + seed[j] + j; /* non linear */
         i++; j++;
         if (i >= NN) { _state[0] = _state[NN-1]; i = 1; }
         if (j >= seed_size) { j = 0; }
      }
      for (k = NN-1; k > 0; k--) {
         _state[i] = (_state[i] ^ ((_state[i-1] ^ (_state[i-1] >> 62)) * 2862933555777941757ULL))
            - i; /* non linear */
         i++;
         if (i >= NN) { _state[0] = _state[NN-1]; i = 1; }
      }
     _state[0] = 1ULL << 63; /* MSB is 1; assuring non-zero initial array */ 
   }
}

/*
 * Copy constructor.
 * Give the copy the same state as the original.
 */
mersenne_twister_64::mersenne_twister_64(const mersenne_twister_64& r)
 : _state(r._state),
   _curr(r._curr)
{ }

/*
 * Destructor.
 */
mersenne_twister_64::~mersenne_twister_64() {
   /* do nothing */
}

/*
 * Generate a pseudorandom unsigned 64-bit integer uniformly distributed 
 * on the closed interval [0, 2^64-1].
 */
unsigned long long mersenne_twister_64::generate() {
   static const unsigned long long mag01[2] = {0ULL, MATRIX_A};
   unsigned long long x;
   /* check if new items need to be generated */
   if (_curr >= NN) {
      /* generate NN items at once */
      unsigned long i;
      for (i = 0; i < (NN-MM); i++) {
         x = (_state[i] & UM) | (_state[i+1] & LM);
         _state[i] = _state[i+MM] ^ (x >> 1) ^ mag01[static_cast<unsigned long>(x & 1ULL)];
      }
      for (; i < (NN-1); i++) {
         x = (_state[i] & UM) | (_state[i+1] & LM);
         _state[i] = _state[i + MM - NN] ^ (x >> 1) ^ mag01[static_cast<unsigned long>(x & 1ULL)];
      }
      x = (_state[NN-1] & UM) | (_state[0] & LM);
      _state[NN-1] = _state[MM-1] ^ (x >> 1) ^ mag01[static_cast<unsigned long>(x & 1ULL)];
      _curr = 0;
   }
   /* retrieve next item */
   x = _state[_curr++];
   x ^= (x >> 29) & 0x5555555555555555ULL;
   x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
   x ^= (x << 37) & 0xFFF7EEE000000000ULL;
   x ^= (x >> 43);
   return x;
}

} /* namespace sources */
} /* namespace random */
} /* namespace math */
