// Copyright (C) 2002 Charless C. Fowlkes <fowlkes@eecs.berkeley.edu>
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

#include <iostream>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <stdarg.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include "exception.hh"
#include "string.hh"
#include "util.hh"
#include "array.hh"
#include "message.hh"

namespace Util {

    void
     registerConfig() {
        Message::registerConfig();
        Random::registerConfig();
    } void

     init() {
        Message::init();
        Random::init();
    }

    FILE *openInputStrm(const char *fileName) throw(Exception) {
        if (strcmp(fileName, "-") == 0) {
            return stdin;
        }
        if (strcmp(fileName, "stdin") == 0) {
            return stdin;
        }
        FILE *strm = fopen(fileName, "r");
        if (strm == NULL) {
            throw
                Exception(String
                          ("Could not open file \"%s\" for reading.\n",
                           fileName));
        }
        return strm;
    }

    FILE *openOutputStrm(const char *fileName) throw(Exception) {
        if (strcmp(fileName, "-") == 0) {
            return stdout;
        }
        if (strcmp(fileName, "stdout") == 0) {
            return stdout;
        }
        if (strcmp(fileName, "stderr") == 0) {
            return stderr;
        }
        FILE *strm = fopen(fileName, "w");
        if (strm == NULL) {
            throw
                Exception(String
                          ("Could not open file \"%s\" for writing.\n",
                           fileName));
        }
        return strm;
    }

    void
     bucketSort(const uint * data, uint nElem, uint maxVal,
                uint * perm, uint * scratch) {
        uint *counts;
        if (scratch == NULL) {
            counts = new uint[maxVal + 1];
        } else {
            counts = scratch;
        }

        // Initialize counts.
        for (uint i = 0; i <= maxVal; i++) {
            counts[i] = 0;
        }

        // Count the values.
        for (uint i = 0; i < nElem; i++) {
            uint val = data[i];
            assert(val <= maxVal);
            counts[val]++;
        }

        // Compute parallel prefix of counts, shifted right one element.
        for (uint sum = 0, i = 0; i <= maxVal; i++) {
            uint tmp = counts[i];
            counts[i] = sum;
            sum += tmp;
        }

        // Compute permutation vector.
        for (uint i = 0; i < nElem; i++) {
            uint val = data[i];
            assert(val <= maxVal);
            assert(counts[val] < nElem);
            perm[counts[val]++] = i;
        }

        // Delete scratch space if we allocated it.
        if (scratch == NULL) {
            delete[]counts;
        }
    }

    static const uint b = 10;

    uint bucketSortCompressScratchSize(uint maxVal) {
        return (maxVal >> (b + 12)) + 1 + (maxVal >> 6) + 1;
    }

    void
     bucketSortCompress(uint * data, uint & nElem, uint maxVal,
                        uint64 * scratch) {
        uint64 *smallBitvec;
        uint64 *largeBitvec;
        uint n = bucketSortCompressScratchSize(maxVal);
        if (scratch == NULL) {
            uint64 *scratch = new uint64[n];
            smallBitvec = &scratch[0];
            largeBitvec = &scratch[(maxVal >> (b + 12)) + 1];
        } else {
            smallBitvec = &scratch[0];
            largeBitvec = &scratch[(maxVal >> (b + 12)) + 1];
        }

        // Initialize bitvector to zero.
        memset(smallBitvec, 0, n * 8);

        // Set bits based on the values.
        for (uint i = 0; i < nElem; i++) {
            uint val = data[i];
            assert(val <= maxVal);
            uint index, offset;
            index = val >> (b + 12);
            offset = (val >> (b + 6)) & 0x3f;
            smallBitvec[index] |= 1 << offset;
            index = val >> 6;
            offset = val & 0x3f;
            largeBitvec[index] |= 1 << offset;
        }

        // Find the set bits.
        uint count = 0;
        for (uint index = 0; index < ((maxVal >> (b + 12)) + 1); index++) {
            uint64 word = smallBitvec[index];
            if (word == 0) {
                continue;
            }
            for (uint offset = 0; offset < 64; offset++) {
                if ((word & (1 << offset)) == 0) {
                    continue;
                }
                uint start = ((index << 6) + offset) << b;
                for (uint index = start; index < start + (1 << b); index++) {
                    uint64 word = largeBitvec[index];
                    if (word == 0) {
                        continue;
                    }
                    for (uint offset = 0; offset < 64; offset++) {
                        if ((word & (1 << offset)) == 0) {
                            continue;
                        }
                        uint val = (index << 6) + offset;
                        data[count++] = val;
                    }
                }
            }
        }
        nElem = count;

        // Delete scratch space if we allocated it.
        if (scratch == NULL) {
            delete[]smallBitvec;
        }
    }

// O(n) implementation.
    static void
     _kOfN_largeK(uint k, uint n, uint * values) {
        assert(k > 0);
        assert(k <= n);
        uint j = 0;
        for (uint i = 0; i < n; i++) {
            double prob = (double) (k - j) / (n - i);
            assert(prob <= 1);
            double x = Util::rand->fp();
            if (x < prob) {
                values[j++] = i;
            }
        }
        assert(j == k);
    }

// O(k*lg(k)) implementation; constant factor is about 2x the constant
// factor for the O(n) implementation.
    static void
     _kOfN_smallK(uint k, uint n, uint * values) {
        assert(k > 0);
        assert(k <= n);
        if (k == 1) {
            values[0] = Util::rand->ui32(0, n - 1);
            return;
        }
        uint leftN = n / 2;
        uint rightN = n - leftN;
        uint leftK = 0;
        uint rightK = 0;
        for (uint i = 0; i < k; i++) {
            uint x = Util::rand->ui32(0, n - i - 1);
            if (x < leftN - leftK) {
                leftK++;
            } else {
                rightK++;
            }
        }
        if (leftK > 0) {
            _kOfN_smallK(leftK, leftN, values);
        }
        if (rightK > 0) {
            _kOfN_smallK(rightK, rightN, values + leftK);
        }
        for (uint i = leftK; i < k; i++) {
            values[i] += leftN;
        }
    }

// Return k randomly selected integers from the interval [0,n), in
// increasing sorted order.
    void
     kOfN(uint k, uint n, uint * values) {
        if (k == 0) {
            return;
        }
        static double log2 = log(2);
        double klogk = k * log(k) / log2;
        if (klogk < n / 2) {
            _kOfN_smallK(k, n, values);
        } else {
            _kOfN_largeK(k, n, values);
        }
    }

    void
     bug(const char *fmt, ...) {
        va_list ap;

        fflush(stdout);

        fprintf(stderr, "\nINTERNAL BUG: ");
        va_start(ap, fmt);
        vfprintf(stderr, fmt, ap);
        va_end(ap);

        abort();
    }

    void
     fatal(const char *fmt, ...) {
        va_list ap;

        fflush(stdout);

        fprintf(stderr, "\nFATAL ERROR: ");
        va_start(ap, fmt);
        vfprintf(stderr, fmt, ap);
        va_end(ap);

        exit(1);
    }

    void
     warning(const char *fmt, ...) {
        va_list ap;

        fflush(stdout);

        fprintf(stderr, "WARNING: ");
        va_start(ap, fmt);
        vfprintf(stderr, fmt, ap);
        va_end(ap);
    }

    void
     system(const char *cmd) {
        std::cerr << "Executing command: " << cmd << std::endl;
        int status =::system(cmd);
        if (status == -1) {
            throw Exception(String("Util::sytem(): system() failed."));
        }
        if (status != 0) {
            throw Exception(String("Util::sytem(%s): status = %d.",
                                   cmd, status));
        }
    }

// Execute the command.
// Argument list must be NULL-terminated!
// Return only if the command exited normally.
    void
     execute(const char *path, ...) {
        va_list ap;

        // count the number of arguments
        int argc = 1;
        va_start(ap, path);
        while (va_arg(ap, const char *) != NULL) {
            argc++;
        }
        va_end(ap);

        // create an array of arguments
        Util::Array1D < const char *>argv(argc + 1);
        va_start(ap, path);
        argv(0) = path;
        String command("%s", argv(0));
        for (int i = 1; i <= argc; i++) {
            argv(i) = va_arg(ap, const char *);
            if (i < argc) {
                command.append(" %s", argv(i));
            }
        }
        va_end(ap);
        assert(argv(argc) == NULL);

        std::cerr << "Executing command: " << command << std::endl;

        // fork and execute command
        int pid = fork();
        if (pid == -1) {
            throw Exception(String("Util::system(): fork() failed."));
        }
        if (pid == 0) {
            // child : exec the command
            (void) execvp(path, (char *const *) argv.data());
            // execvp() should never return
            throw Exception(String("Util::system(): execvp(%s,%s) failed.",
                                   path, command.text()));
        } else {
            // parent : wait for child
            int status;
            int rval = waitpid(pid, &status, 0);
            // make sure waitpid() call succeeded
            if (rval != pid) {
                throw
                    Exception(String("Util::sytem(): waitpid() failed."));
            }
            // make sure child exited normally, under its own
            // control
            if (!WIFEXITED(status)) {
                throw
                    Exception(String
                              ("Util::sytem(): command '%s' exited abnormally.",
                               command.text()));
            }
            // make sure child's exit status was 0
            if (WEXITSTATUS(status) != 0) {
                throw
                    Exception(String
                              ("Util::sytem(): command '%s' exited with status=%d.",
                               command.text(), WEXITSTATUS(status)));
            }
            // command succeeded
            return;
        }
    }

// execvp the program, and return the ends of the pipes corresponding
// to the program's stdin, stdout, stderr.
    int
     systemPipe(const char *file, char *const argv[],
                FILE * &in, FILE * &out) {
        // create the pipes
        int inPipe[2];
        int outPipe[2];
        if (pipe(inPipe) == -1) {
            throw Exception("systemPipe: pipe failed");
        }
        if (pipe(outPipe) == -1) {
            throw Exception("systemPipe: pipe failed");
        }
        // get process group id of the parent
        pid_t pgid = getpgrp();
        if (pgid == -1) {
            throw Exception("systemPipe: getpgid failed");
        }
        // fork
        pid_t pid = fork();
        if (pid == -1) {
            throw Exception("systemPipe: fork failed");
        }

        if (pid == 0) {         // CHILD
            // set process group id to that of the parent
            if (setpgid(0, pgid) == -1) {
                std::cerr << "systemPipe: setpgid failed" << std::endl;
                exit(1);
            }
            // bind the pipe ends to stdin and stdout
            if (dup2(inPipe[0], 0) == -1) {
                std::cerr << "systemPipe: dup2 failed" << std::endl;
                exit(1);
            }
            if (dup2(outPipe[1], 1) == -1) {
                std::cerr << "systemPipe: dup2 failed" << std::endl;
                exit(1);
            }
            // close unused file descriptors
            if (close(inPipe[1]) == -1) {
                std::cerr << "systemPipe: close failed" << std::endl;
                exit(1);
            }
            if (close(outPipe[0]) == -1) {
                std::cerr << "systemPipe: close failed" << std::endl;
                exit(1);
            }
            if (close(inPipe[0]) == -1) {
                std::cerr << "systemPipe: close failed" << std::endl;
                exit(1);
            }
            if (close(outPipe[1]) == -1) {
                std::cerr << "systemPipe: close failed" << std::endl;
                exit(1);
            }
            // exec
            if (execvp(file, argv) == -1) {
                std::cerr << "systemPipe: execvp failed" << std::endl;
                exit(1);
            }
        }
        // close unused file descriptors
        if (close(inPipe[0]) == -1) {
            throw Exception("systemPipe: close failed");
        }
        if (close(outPipe[1]) == -1) {
            throw Exception("systemPipe: close failed");
        }
        // recast pipe ends as streams
        in = fdopen(inPipe[1], "w");
        out = fdopen(outPipe[0], "r");
        if (in == NULL || out == NULL) {
            throw Exception("systemPipe: fdopen failed");
        }
        // return pid of child
        return pid;
    }

    double mod (double x, double y) {
      return (y==0) ? x : (x - y*floor(x/y));
    }

}                              // namespace Util
