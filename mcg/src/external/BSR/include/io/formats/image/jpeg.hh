/*
 * Jpeg file utilities.
 */
#ifndef IO__FORMATS__IMAGE__JPEG_HH
#define IO__FORMATS__IMAGE__JPEG_HH

#include "lang/pointers/auto_ptr.hh"
#include "lang/string.hh"
#include "math/matrices/matrix.hh"

namespace io {
namespace formats {
namespace image {
/*
 * Imports.
 */
using lang::pointers::auto_ptr;
using lang::string;
using math::matrices::matrix;

/*
 * Jpeg file utilities.
 */
class jpeg {
public:
   /*
    * Check whether the given jpeg file contains a grayscale or color image.
    */
   static bool is_grayscale(const string<>& /* filename */);

   /*
    * Read a grayscale image from a jpeg file.
    * The returned image has values in [0,1].
    * Throw an exception if the file contains a color image.
    */
   static void read(
      const string<>&,        /* filename */
      auto_ptr< matrix<> >&   /* returned grayscale image */
   );

   /*
    * Read an RGB color image from a jpeg file.
    * The returned red, green, and blue channels have values in [0,1].
    * Throw an exception if the file contains a grayscale image.
    */
   static void read(
      const string<>&,        /* filename */
      auto_ptr< matrix<> >&,  /* returned red channel */
      auto_ptr< matrix<> >&,  /* returned green channel */
      auto_ptr< matrix<> >&   /* returned blue channel */
   );
      
   /*
    * Write a grayscale image to a jpeg file.
    * The image should have values in [0,1].
    * Optionally specify a compression quality in [0,1].
    */
   static void write(
      const string<>&,  /* filename */
      const matrix<>&,  /* grayscale image */
      double = 1.0      /* quality (default of one is best image quality) */
   );

   /* 
    * Write an RGB color image to a jpeg file.
    * The red, green, and blue channels should have values in [0,1].
    * Optionally specify a compression quality in [0,1].
    */
   static void write(
      const string<>&,  /* filename */
      const matrix<>&,  /* red channel */
      const matrix<>&,  /* green channel */
      const matrix<>&,  /* blue channel */
      double = 1.0      /* quality (default of one is best image quality) */
   );
};
 
} /* namespace image */
} /* namespace formats */
} /* namespace io */

#endif
