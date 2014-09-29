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


#include <sys/time.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "random.hh"
#include "message.hh"
#include "string.hh"
#include "exception.hh"
#include "configure.hh"

namespace Util {

    Random *rand = NULL;

    void Random::registerConfig() {
        static bool called = false;
        if (called) {
            return;
        }
        Configure::registerInt("Random::seed", 0,
                               "Seed for the random number generator.  If zero, then the seed is "
                               "taken from the system clock.");
        called = true;
    }

    void Random::init() {
        static bool called = false;
        if (called) {
            return;
        }
        uint64 seed = Configure::getInt("Random::seed");
        if (seed == 0) {
            struct timeval t;
            gettimeofday(&t, NULL);
            uint64 a = (t.tv_usec >> 3) & 0xffff;
            uint64 b = t.tv_sec & 0xffff;
            uint64 c = (t.tv_sec >> 16) & 0xffff;
            seed = a | (b << 16) | (c << 32);
        }
        seed &= 0xffffffffffffull;
        assert(Util::rand == NULL);
        Util::rand = new Random(seed);
        Message::debug(String("random seed: 0x%llx", seed),1);
        called = true;
    }

    Random *Random::spawn() {
        return new Random(this);
    }

    void
     Random::reset() {
        _init(_seed);
    }

    Random::Random(uint64 seed) {
        _init(seed);
    }

    Random::Random(Random * that) {
        uint64 a = that->ui32();
        uint64 b = that->ui32();
        uint64 seed = (a << 32) | b;
        _init(seed);
    }

    void
     Random::_init(uint64 seed) {
        _seed = seed & 0xffffffffffffull;
        _xsubi[0] = (seed >> 0) & 0xffff;
        _xsubi[1] = (seed >> 16) & 0xffff;
        _xsubi[2] = (seed >> 32) & 0xffff;
    }

}                              // namespace Util
