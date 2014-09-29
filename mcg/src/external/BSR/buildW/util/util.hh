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

//#ifndef UITL_HH
#ifndef UTIL_HH // 05/2014 - Modified by Jordi Pont-Tuset <jordi.pont@upc.edu>
#define UTIL_HH

#include <stdio.h>
#include <ctype.h>
#include "types.hh"
#include "random.hh"
#include "exception.hh"
#include "sort.hh"

namespace Util {

    void registerConfig();
    void init();

    // Creates permutation vector, where the sorted data is given
    // by data[perm[i]] for i=[0,nElem).
    // Requires O(maxVal) temporary space.
    // Takes O(nElem+maxVal) time.
    // Scratch array must have at least maxVal+1 elements, if provided.
    void bucketSort(const uint * data, uint nElem, uint maxVal,
                    uint * perm, uint * scratch = NULL);

    // Return k unique values chosen from [0..n), in sorted order.
    // O(min{n,k*lg(k)}) time.
    void kOfN(uint k, uint n, uint * values);

    // WARNING: this routine is not debugged at all!
    // Overwrites data/nElem with sorted data, duplicates removed.
    // If scratch space is provided, it must be at least as big as
    // the value returned by bucketSortCompressScratchSize().
    // Requires O(maxVal) space and time.
    // Uses a 2-level bit vector; thus, this algorithm should be
    // fast when nElem << maxVal, and the values are clustered.
    uint bucketSortCompressScratchSize(uint maxVal);
    void bucketSortCompress(uint * data, uint & nElem, uint maxVal,
                            uint64 * scratch = NULL);

    // Safe file opening.
    // Filename may be {-,stdin,stdout,stderr} where appropriate.
    FILE *openInputStrm(const char *fileName) throw(Exception);
    FILE *openOutputStrm(const char *fileName) throw(Exception);

    // Common mini functions.
#undef abs
#undef min
#undef max

    #define EPSILON 0.00000000001
     
    //wrap values in [-max,2*max) into [0,max)
    //dangerous when instantiated as float/double!
    template < class T > inline T wrap(const T val, const T max)
    {
        T rval = val;
        if (rval < 0)
          {
              rval += max;
          }
        if (rval >= max)
          {
              rval -= max;
          }
        assert (rval >= 0 && rval < max);
        return rval;
    }
     
    double mod (double x, double y);

    template < class T > T abs(T a) {
        return a < 0 ? -a : a;
    } 
    template < class T > T max(T a, T b) {
        return a > b ? a : b;
    }
    template < class T > T min(T a, T b) {
        return a < b ? a : b;
    }
    template < class T > T max(T a, T b, T c) {
        return max(a, max(b, c));
    }
    template < class T > T min(T a, T b, T c) {
        return min(a, min(b, c));
    }
    template < class T > T max(T a, T b, T c, T d) {
        return max(a, max(b, c, d));
    }
    template < class T > T min(T a, T b, T c, T d) {
        return min(a, min(b, c, d));
    }
    template < class T > T minmax(T a, T b, T c) {
        return max(a, min(b, c));
    }

    template < class T > void swap(T & a, T & b) {
        T tmp(a);
        a = b;
        b = tmp;
    }

    inline double gaussian(double sigma, double x) {
        return 1.0 / (sqrt(2 * M_PI) * sigma) * exp(-(x * x) / (2 * sigma * sigma));
    }

    // Arithmetic/Geometric/Harmonic mean.
    template < class T > double aMean(T * a, int n) {
        if (n == 0) {
            return 0;
        }
        double sum = 0;
        for (int i = 0; i < n; i++) {
            sum += a[i];
        }
        return sum / n;
    }
    template < class T > double gMean(T * a, int n) {
        if (n == 0) {
            return 0;
        }
        double sum = 0;
        for (int i = 0; i < n; i++) {
            sum += log(a[i]);
        }
        return pow(M_E, sum / n);
    }
    template < class T > double hMean(T * a, int n) {
        if (n == 0) {
            return 0;
        }
        double sum = 0;
        for (int i = 0; i < n; i++) {
            sum += 1.0 / a[i];
        }
        return n / sum;
    }

    // sample variance
    template < class T > double variance(T * a, int n) {
        if (n <= 1) {
            return 0;
        }
        const double mu = aMean(a, n);
        double var = 0;
        for (int i = 0; i < n; i++) {
            const double v = a[i] - mu;
            var += v * v;
        }
        var /= n - 1;
        return var;
    }

    // Randomize the elements of an array.
    template < class T > void randomizeArray(uint n, T * a) {
        // For i=[0..n), move random node in [i..n) into slot i.
        for (uint i = 0; i < n; i++) {
            uint j = Util::rand->ui32(i, n - 1);
            swap(a[i], a[j]);
        }
    }

    // Compute the difference between two angles.
    // t1,t2 in [0,PI)
    // diff in [0,PI/2]
    inline double
     thetaDiff(double t1, double t2) {
        assert(finite(t1));
        assert(finite(t2));
        assert(t1 >= 0 && t1 < M_PI + 1e-5);
        assert(t2 >= 0 && t2 < M_PI + 1e-5);
        double d = fabs(t1 - t2);
        if (d > M_PI_2) {
            d = M_PI - d;
        }
        if (d > M_PI_2) {
            d = M_PI_2;
        }
        if (d < 0) {
            d = 0;
        }
        assert(finite(d));
        assert(d >= 0 && d <= M_PI_2);
        return d;
    }

    // Execute the command.  
    // Argument list for execute() must be NULL-terminated!
    // Return only if the command exited normally.
    void system(const char *cmd);
    void execute(const char *path, ...);

    // execvp the program, and return the ends of the pipes corresponding
    // to the program's stdin, stdout, stderr.
    // Return the PID of the child.
    int systemPipe(const char *file, char *const argv[],
                   FILE * &in, FILE * &out);

    // Call when an internal bug is detected.
    // Immediately stops program via abort().
    void bug(const char *fmt, ...);

    // Call when a fatal error is detected.
    // Immediately stops program via exit().
    void fatal(const char *fmt, ...);

    // Call to report a run-time warning.
    // Does not stop program.
    void warning(const char *fmt, ...);

}                              // namespace Util

#endif                          // __Util_h__
