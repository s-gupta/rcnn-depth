// Copyright (C) 2002 David R. Martin <dmartin@eecs.berkeley.edu>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of the
// License, or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
// 02111-1307, USA, or see http://www.gnu.org/copyleft/gpl.html.

#ifndef RANDOM_HH
#define RANDOM_HH

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "types.hh"

namespace Util 
{

    // All random numbers are generated from a single seed.  This is true
    // even when private random streams (seperate from the global
    // Util::rand stream) are spawned from existing streams, since the new
    // streams are seeded automatically from the parent's random stream.
    // Any random stream can be reset so that a sequence of random values
    // can be replayed.

    // If seed==0, then the seed is generated from the system clock.  The
    // init() method will always print the actual seed to stderr.

    class Random 
    {
      public:

        static void registerConfig();
        static void init();

        // Spawn off a separate random stream.  The stream is seeded
        // with the parent's random stream.
         Random(Random * that);
        Random *spawn();

        // Restore initial seed so we can replay a random sequence.
        void reset();

        // double in [0..1) or [a..b)
        inline double fp();
        inline double fp(double a, double b);

        // int32 in [INT32_MIN..INT32_MAX] or [a..b]
        inline int32 i32();
        inline int32 i32(int32 a, int32 b);

        // uint32 in [0..UINT32_MAX] or [a..b]
        inline uint32 ui32();
        inline uint32 ui32(uint32 a, uint32 b);

      protected:

         Random(uint64 seed);
        void _init(uint64 seed);

        // The original seed for this random stream.
        uint64 _seed;

        // The current state for this random stream.
        uint16 _xsubi[3];

    };

    extern Random *rand;

    inline uint32 Random::ui32() {
        // uniform uint32 in [0..2^32)
        return ui32(0, UINT32_MAX);
    } 
    
    inline uint32 Random::ui32(uint32 a, uint32 b) {
        // uniform uint32 in [a..b]
        assert(a <= b);
        double x = fp();
        return (uint32) floor(x * ((double) b - (double) a + 1) + a);
    }

    inline int32 Random::i32() {
        // uniform int32 in [-2^31..2^31)
        return i32(INT32_MIN, INT32_MAX);
    }

    inline int32 Random::i32(int32 a, int32 b) {
        // uniform int32 in [a..b]
        assert(a <= b);
        double x = fp();
        return (int32) floor(x * ((double) b - (double) a + 1) + a);
    }

    inline double Random::fp() {
        // uniform double in [0..1)
        return erand48(_xsubi);
    }

    inline double Random::fp(double a, double b) {
        // uniform double in [a..b)
        assert(a < b);
        return erand48(_xsubi) * (b - a) + a;
    }

}                              // namespace Util

#endif                          // __Random_h__
