/*
 * Pseudorandom source.
 * A pseudorandom source is a generator that serves as a source of random 
 * numbers for all of the built-in numeric types.
 */
#ifndef MATH__RANDOM__SOURCES__RAND_SOURCE_HH
#define MATH__RANDOM__SOURCES__RAND_SOURCE_HH

namespace math {
namespace random {
namespace sources {

/*
 * Abstract base class for pseudorandom source.
 */
class rand_source {
public:
   /*
    * Destructor.
    */
   virtual ~rand_source() = 0;

   /*
    * Generate a pseudorandom boolean uniformly distributed on {true, false}.
    */
   virtual bool gen_uniform_bool() = 0;
   
   /*
    * Generate a pseudorandom integer uniformly distributed on the closed
    * interval [0, 2^(N-1) - 1], where N is the number of bits of precision in
    * the corresponding unsigned datatype.
    */
   virtual char      gen_uniform_char();
   virtual short     gen_uniform_short();
   virtual int       gen_uniform_int();
   virtual long      gen_uniform_long();
   virtual long long gen_uniform_long_long();

   /*
    * Generate a pseudorandom integer uniformly distributed on the closed
    * interval [0, 2^N - 1], where N is the number of bits of precision in
    * the corresponding datatype.
    */
   virtual unsigned char      gen_uniform_unsigned_char() = 0;
   virtual unsigned short     gen_uniform_unsigned_short() = 0;
   virtual unsigned int       gen_uniform_unsigned_int() = 0;
   virtual unsigned long      gen_uniform_unsigned_long() = 0;
   virtual unsigned long long gen_uniform_unsigned_long_long() = 0;

   /*
    * Generate a pseudorandom real number uniformly distributed on the 
    * closed interval [0,1].
    */
   virtual float       gen_uniform_closed_float();
   virtual double      gen_uniform_closed_double();
   virtual long double gen_uniform_closed_long_double();
   
   /*
    * Generate a pseudorandom real number uniformly distributed on the
    * half-open interval [0,1).
    */
   virtual float       gen_uniform_half_open_float();
   virtual double      gen_uniform_half_open_double();
   virtual long double gen_uniform_half_open_long_double();

   /*
    * Generate a pseudorandom real number uniformly distributed on the
    * open interval (0,1).
    */
   virtual float       gen_uniform_open_float();
   virtual double      gen_uniform_open_double();
   virtual long double gen_uniform_open_long_double();

   /*
    * Generate a pseudorandom real number normally distributed with 
    * mean zero and unit variance.
    */
   float       gen_normal_float();
   double      gen_normal_double();
   long double gen_normal_long_double();

   /*
    * Generate a pseudorandom real number normally distributed with the
    * specified mean and variance.
    */
   float gen_normal_float(
      float,         /* mean */
      float          /* variance */
   );
   
   double gen_normal_double(
      double,        /* mean */
      double         /* variance */
   );
   
   long double gen_normal_long_double(
      long double,   /* mean */
      long double    /* variance */
   );
};

} /* namespace sources */
} /* namespace random */
} /* namespace math */

#endif
