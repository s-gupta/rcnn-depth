/*
 * Arrays.
 */
#include "lang/array.hh"

namespace lang {

/*
 * Array initialization functions for built-in types.
 * These functions guarantee that arrays are always zero-initialized.
 */
void array_initialize(unsigned long size, bool* data) {
   for (unsigned long n = 0; n < size; n++)
      data[n] = false;
}
 
void array_initialize(unsigned long size, char* data) {
   for (unsigned long n = 0; n < size; n++)
      data[n] = 0;
}

void array_initialize(unsigned long size, unsigned char* data) {
   for (unsigned long n = 0; n < size; n++)
      data[n] = 0;
}

void array_initialize(unsigned long size, short* data) {
   for (unsigned long n = 0; n < size; n++)
      data[n] = 0;
}

void array_initialize(unsigned long size, unsigned short* data) {
   for (unsigned long n = 0; n < size; n++)
      data[n] = 0;
}

void array_initialize(unsigned long size, int* data) {
   for (unsigned long n = 0; n < size; n++)
      data[n] = 0;
}

void array_initialize(unsigned long size, unsigned int* data) {
   for (unsigned long n = 0; n < size; n++)
      data[n] = 0;
}

void array_initialize(unsigned long size, long* data) {
   for (unsigned long n = 0; n < size; n++)
      data[n] = 0;
}

void array_initialize(unsigned long size, unsigned long* data) {
   for (unsigned long n = 0; n < size; n++)
      data[n] = 0;
}

void array_initialize(unsigned long size, long long* data) {
   for (unsigned long n = 0; n < size; n++)
      data[n] = 0;
}

void array_initialize(unsigned long size, unsigned long long* data) {
   for (unsigned long n = 0; n < size; n++)
      data[n] = 0;
}

void array_initialize(unsigned long size, float* data) {
   for (unsigned long n = 0; n < size; n++)
      data[n] = 0;
}

void array_initialize(unsigned long size, double* data) {
   for (unsigned long n = 0; n < size; n++)
      data[n] = 0;
}

void array_initialize(unsigned long size, long double* data) {
   for (unsigned long n = 0; n < size; n++)
      data[n] = 0;
}

} /* namespace lang */
