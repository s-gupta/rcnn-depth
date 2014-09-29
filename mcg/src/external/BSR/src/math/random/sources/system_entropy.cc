/*
 * Pseudorandom source that uses system entropy to generate random numbers.
 */
#include "concurrent/threads/synchronization/locks/auto_lock.hh"
#include "concurrent/threads/synchronization/mutex.hh"
#include "io/streams/ifstream.hh"
#include "lang/types/type_sizes.hh"
#include "math/random/sources/system_entropy.hh"

namespace math {
namespace random {
namespace sources {
/*
 * Imports.
 */
using concurrent::threads::synchronization::locks::auto_lock;
using concurrent::threads::synchronization::mutex;
 
/*
 * Anonymous namespace for random device.
 */
namespace {
using concurrent::threads::synchronization::mutex;
using io::streams::ifstream;
/*
 * Initialize mutex guarding random device to unlocked state.
 */
static mutex rand_dev_mutex;

/*
 * Return the file stream for the random device.
 */
ifstream& rand_dev_stream() {
   static ifstream* f = new ifstream(
      "/dev/urandom", 
      ifstream::in | ifstream::binary
   );
   return *f;
}
} /* anonymous namespace for random device */

/*
 * Constructor.
 */
system_entropy::system_entropy() {
   /* do nothing */
}

/*
 * Copy constructor.
 */
system_entropy::system_entropy(const system_entropy&) {
   /* do nothing */
}

/*
 * Destructor.
 */
system_entropy::~system_entropy() {
   /* do nothing */
}

/*
 * Generate a pseudorandom boolean uniformly distributed on {true, false}.
 */
bool system_entropy::gen_uniform_bool() {
   return static_cast<bool>(
      this->gen_uniform_unsigned_char() >> (UCHAR_BITS - 1)
   );
}

/*
 * Generate a pseudorandom integer uniformly distributed on the closed
 * interval [0, 2^N - 1], where N is the number of bits of precision in
 * the corresponding datatype.
 */
unsigned char system_entropy::gen_uniform_unsigned_char() {
   unsigned char x;
   auto_lock<mutex> xlock(rand_dev_mutex);
   rand_dev_stream().read(
      static_cast<char*>(static_cast<void*>(&x)),
      sizeof(unsigned char)
   );
   return x;
}

unsigned short system_entropy::gen_uniform_unsigned_short() {
   unsigned short x;
   auto_lock<mutex> xlock(rand_dev_mutex);
   rand_dev_stream().read(
      static_cast<char*>(static_cast<void*>(&x)),
      sizeof(unsigned short)
   );
   return x;
}

unsigned int system_entropy::gen_uniform_unsigned_int() {
   unsigned int x;
   auto_lock<mutex> xlock(rand_dev_mutex);
   rand_dev_stream().read(
      static_cast<char*>(static_cast<void*>(&x)),
      sizeof(unsigned int)
   );
   return x;
}

unsigned long system_entropy::gen_uniform_unsigned_long() {
   unsigned long x;
   auto_lock<mutex> xlock(rand_dev_mutex);
   rand_dev_stream().read(
      static_cast<char*>(static_cast<void*>(&x)),
      sizeof(unsigned long)
   );
   return x;
}

unsigned long long system_entropy::gen_uniform_unsigned_long_long() {
   unsigned long long x;
   auto_lock<mutex> xlock(rand_dev_mutex);
   rand_dev_stream().read(
      static_cast<char*>(static_cast<void*>(&x)),
      sizeof(unsigned long long)
   );
   return x;
}

} /* namespace sources */
} /* namespace random */
} /* namespace math */
