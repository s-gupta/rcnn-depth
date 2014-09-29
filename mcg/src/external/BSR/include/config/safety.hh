/*
 * Safety settings.
 */
#ifndef CONFIG__SAFETY_HH
#define CONFIG__SAFETY_HH

/*
 * Check smart pointer dereference?
 * If true, smart pointers should throw exceptions on null dereference.
 */
#define CONFIG__SAFETY__CHECK_DEREFERENCE (true)

/*
 * Check bounds when accessing arrays?
 * If true, arrays should throw exceptions on access to an out of bounds index.
 */
#define CONFIG__SAFETY__CHECK_BOUNDS (true)

#endif
