/*
 * Ranges of built-in data types.
 */
#ifndef LANG__TYPES__TYPE_RANGES_HH
#define LANG__TYPES__TYPE_RANGES_HH

#include <cfloat>
#include <climits>

/*
 * Range of built-in integer types.
 *
 * These are system-specific values.  Any definitions commented out are 
 * replaced by system-specific values defined in the included headers (and 
 * hence commented out values are for illustrative purposes only).
 */
//#define CHAR_MIN                           -128
//#define SHRT_MIN                         -32768
//#define INT_MIN                     -2147483648
//#define LONG_MIN                    -2147483648L  (32-bit systems)
//#define LONG_MIN           -9223372036854775808L  (64-bit systems)
//#define LONG_LONG_MIN      -9223372036854775808LL

//#define CHAR_MAX                            127
//#define SHRT_MAX                          32767
//#define INT_MAX                      2147483647
//#define LONG_MAX                     2147483647L  (32-bit systems)
//#define LONG_MAX            9223372036854775807L  (64-bit systems)
//#define LONG_LONG_MAX       9223372036854775807LL

#define UCHAR_MIN                               0
#define USHRT_MIN                               0
#define UINT_MIN                                0
#define ULONG_MIN                               0
#define ULONG_LONG_MIN                          0

//#define UCHAR_MAX                           255
//#define USHRT_MAX                         65535
//#define UINT_MAX                     4294967295
//#define ULONG_MAX                    4294967295UL (32-bit systems)
//#define ULONG_MAX          18446744073709551615UL (64-bit systems)
#ifdef __APPLE__
#define ULONG_LONG_MAX     18446744073709551615ULL
#endif

/*
 * Range and precision of built-in floating point types.
 *
 * These are system-specific values.  Any definitions commented out are 
 * replaced by system-specific values defined in the included headers (and 
 * hence commented out values are for illustrative purposes only).
 */
//#define FLT_RADIX                             2
//#define FLT_MANT_DIG                         24
//#define FLT_DIG                               6
//#define FLT_MIN_EXP                        -125
//#define FLT_MIN_10_EXP                      -37
//#define FLT_MAX_EXP                         128
//#define FLT_MAX_10_EXP                      +38
//#define FLT_MIN                  1.17549435E-38F
//#define FLT_MAX                  3.40282347E+38F
//#define FLT_EPSILON              1.19209290E-07F

//#define DBL_MANT_DIG                         53
//#define DBL_DIG                              15
//#define DBL_MIN_EXP                       -1021
//#define DBL_MIN_10_EXP                     -307
//#define DBL_MAX_EXP                        1024
//#define DBL_MAX_10_EXP                      308
//#define DBL_MIN         2.2250738585072014E-308
//#define DBL_MAX         1.7976931348623157E+308
//#define DBL_EPSILON     2.2204460492503131E-016

//#define LDBL_MANT_DIG                        64
//#define LDBL_DIG                             18
//#define LDBL_MIN_EXP                     -16381
//#define LDBL_MIN_10_EXP                   -4931
//#define LDBL_MAX_EXP                      16384
//#define LDBL_MAX_10_EXP                    4932
//#define LDBL_MIN       3.3621031431120935E-4932L
//#define LDBL_MAX       1.1897314953572318E+4932L
//#define LDBL_EPSILON     1.0842021724855044E-19L

#endif
