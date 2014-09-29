/*
 * Default pseudorandom source.
 *
 * The current implementation uses the 64-bit mersenne twister.
 */
#include "lang/array.hh"
#include "lang/pointers/auto_ptr.hh"
#include "lang/typecasts/dynamic_typecast.hh"
#include "math/random/sources/rand_source_default.hh"
#include "math/random/sources/mersenne_twister_64.hh"

namespace math {
namespace random {
namespace sources {
/*
 * Impotrs.
 */
using lang::array;
using lang::pointers::auto_ptr;
using lang::typecasts::dynamic_typecast;

/*
 * Constructor.
 * Create a randomly initialized pseudorandom source of the default type.
 */
rand_source_default::rand_source_default()
 : _r(new mersenne_twister_64(array<unsigned long long>()))
{ }

/*
 * Copy constructor.
 */
rand_source_default::rand_source_default(const rand_source_default& r)
 : _r(new mersenne_twister_64(dynamic_typecast<mersenne_twister_64&>(*(r._r))))
{ }

/*
 * Destructor.
 */
rand_source_default::~rand_source_default() {
   /* do nothing */
}

/*
 * Generate a pseudorandom boolean uniformly distributed on {true, false}.
 */
bool rand_source_default::gen_uniform_bool() {
   return _r->gen_uniform_bool();
}

/*
 * Generate a pseudorandom integer uniformly distributed on the closed
 * interval [0, 2^(N-1) - 1], where N is the number of bits of precision in
 * the corresponding unsigned datatype.
 */
char rand_source_default::gen_uniform_char() {
   return _r->gen_uniform_char();
}

short rand_source_default::gen_uniform_short() {
   return _r->gen_uniform_short();
}

int rand_source_default::gen_uniform_int() {
   return _r->gen_uniform_int();
}

long rand_source_default::gen_uniform_long() {
   return _r->gen_uniform_long();
}

long long rand_source_default::gen_uniform_long_long() {
   return _r->gen_uniform_long_long();
}

/*
 * Generate a pseudorandom integer uniformly distributed on the closed
 * interval [0, 2^N - 1], where N is the number of bits of precision in
 * the corresponding datatype.
 */
unsigned char rand_source_default::gen_uniform_unsigned_char() {
   return _r->gen_uniform_unsigned_char();
}

unsigned short rand_source_default::gen_uniform_unsigned_short() {
   return _r->gen_uniform_unsigned_short();
}

unsigned int rand_source_default::gen_uniform_unsigned_int() {
   return _r->gen_uniform_unsigned_int();
}

unsigned long rand_source_default::gen_uniform_unsigned_long() {
   return _r->gen_uniform_unsigned_long();
}

unsigned long long rand_source_default::gen_uniform_unsigned_long_long() {
   return _r->gen_uniform_unsigned_long_long();
}

/*
 * Generate a pseudorandom real number uniformly distributed on the 
 * closed interval [0,1].
 */
float rand_source_default::gen_uniform_closed_float() {
   return _r->gen_uniform_closed_float();
}

double rand_source_default::gen_uniform_closed_double() {
   return _r->gen_uniform_closed_double();
}

long double rand_source_default::gen_uniform_closed_long_double() {
   return _r->gen_uniform_closed_long_double();
}

/*
 * Generate a pseudorandom real number uniformly distributed on the
 * half-open interval [0,1).
 */
float rand_source_default::gen_uniform_half_open_float() {
   return _r->gen_uniform_half_open_float();
}

double rand_source_default::gen_uniform_half_open_double() {
   return _r->gen_uniform_half_open_double();
}

long double rand_source_default::gen_uniform_half_open_long_double() {
   return _r->gen_uniform_half_open_long_double();
}

/*
 * Generate a pseudorandom real number uniformly distributed on the
 * open interval (0,1).
 */
float rand_source_default::gen_uniform_open_float() {
   return _r->gen_uniform_open_float();
}

double rand_source_default::gen_uniform_open_double() {
   return _r->gen_uniform_open_double();
}

long double rand_source_default::gen_uniform_open_long_double() {
   return _r->gen_uniform_open_long_double();
}

} /* namespace sources */
} /* namespace random */
} /* namespace math */
