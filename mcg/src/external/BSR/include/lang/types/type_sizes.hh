/*
 * Sizes of built-in types.
 */
#ifndef LANG__TYPES__TYPE_SIZES_HH
#define LANG__TYPES__TYPE_SIZES_HH

#include <climits>

/*
 * Number of bits per byte.
 */
#define BITS_PER_BYTE      (CHAR_BIT)

/*
 * Sizes of built-in types (in bytes).
 */
#define CHAR_BYTES         (sizeof(char))
#define SHRT_BYTES         (sizeof(short))
#define INT_BYTES          (sizeof(int))
#define LONG_BYTES         (sizeof(long))
#define LONG_LONG_BYTES    (sizeof(long long))

#define UCHAR_BYTES        (sizeof(unsigned char))
#define USHRT_BYTES        (sizeof(unsigned short))
#define UINT_BYTES         (sizeof(unsigned int))
#define ULONG_BYTES        (sizeof(unsigned long))
#define ULONG_LONG_BYTES   (sizeof(unsigned long long))

#define FLT_BYTES          (sizeof(float))
#define DBL_BYTES          (sizeof(double))
#define LDBL_BYTES         (sizeof(long double))

#define PTR_BYTES          (sizeof(void*))

/*
 * Sizes of built-in types (in bits).
 */
#define CHAR_BITS          (BITS_PER_BYTE * CHAR_BYTES)
#define SHRT_BITS          (BITS_PER_BYTE * SHRT_BYTES)
#define INT_BITS           (BITS_PER_BYTE * INT_BYTES)
#define LONG_BITS          (BITS_PER_BYTE * LONG_BYTES)
#define LONG_LONG_BITS     (BITS_PER_BYTE * LONG_LONG_BYTES)

#define UCHAR_BITS         (BITS_PER_BYTE * UCHAR_BYTES)
#define USHRT_BITS         (BITS_PER_BYTE * USHRT_BYTES)
#define UINT_BITS          (BITS_PER_BYTE * UINT_BYTES)
#define ULONG_BITS         (BITS_PER_BYTE * ULONG_BYTES)
#define ULONG_LONG_BITS    (BITS_PER_BYTE * ULONG_LONG_BYTES)

#define FLT_BITS           (BITS_PER_BYTE * FLT_BYTES)
#define DBL_BITS           (BITS_PER_BYTE * DBL_BYTES)
#define LDBL_BITS          (BITS_PER_BYTE * LDBL_BYTES)

#define PTR_BITS           (BITS_PER_BYTE * PTR_BYTES)

#endif
