/*
 * Png file utilities.
 */
#ifndef IO__FORMATS__IMAGE__PNG_HH
#define IO__FORMATS__IMAGE__PNG_HH

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
 * Png file utilities.
 */
class png {
public:
   /*
    * Check whether the given png file contains a grayscale or color image.
    */
   static bool is_grayscale(const string<>& /* filename */);

   /*
    * Check whether the given png file contains an alpha (transparency) channel.
    */
   static bool has_alpha_channel(const string<>& /* filename */);

   /*
    * Read a grayscale image from a png file (discard any alpha channel).
    * The returned image has values in [0,1].
    * Throw an exception if the file contains a color image.
    */
   static void read(
      const string<>&,        /* filename */
      auto_ptr< matrix<> >&   /* returned grayscale image */
   );

   /*
    * Read a grayscale image with alpha channel from a png file.
    * The returned grayscale and alpha channels have values in [0,1].
    * Throw an exception if the file contains a color image.
    *
    * Note that if the file does not define an alpha channel, the
    * returned alpha channel is set to one (completely opaque).
    */
   static void read(
      const string<>&,        /* filename */
      auto_ptr< matrix<> >&,  /* returned grayscale image */
      auto_ptr< matrix<> >&   /* returned alpha channel */
   );

   /*
    * Read an RGB color image from a png file (discard any alpha channel).
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
    * Read an RGBA (RGB with alpha) color image from a png file.
    * The returned red, green, blue, and alpha channels have values in [0,1].
    * Throw an exception if the file contains a grayscale image.
    *
    * Note that if the file does not define an alpha channel, the
    * returned alpha channel is set to one (completely opaque).
    */
   static void read(
      const string<>&,        /* filename */
      auto_ptr< matrix<> >&,  /* returned red channel */
      auto_ptr< matrix<> >&,  /* returned green channel */
      auto_ptr< matrix<> >&,  /* returned blue channel */
      auto_ptr< matrix<> >&   /* returned alpha channel */
   );

   /*
    * Write a grayscale image to a png file.
    * The image should have values in [0,1].
    */
   static void write(
      const string<>&,     /* filename */
      const matrix<>&,     /* grayscale image */
      unsigned long = 16   /* bit depth (8 or 16) */
   );
   
   /*
    * Write a grayscale image with alpha channel to a png file.
    * The grayscale and alpha channels should have values in [0,1].
    */
   static void write(
      const string<>&,     /* filename */
      const matrix<>&,     /* grayscale image */
      const matrix<>&,     /* alpha channel */
      unsigned long = 16   /* bit depth (8 or 16) */
   );

   /* 
    * Write an RGB color image to a png file.
    * The red, green, and blue channels should have values in [0,1].
    */
   static void write(
      const string<>&,     /* filename */
      const matrix<>&,     /* red channel */
      const matrix<>&,     /* green channel */
      const matrix<>&,     /* blue channel */
      unsigned long = 16   /* bit depth (8 or 16) */
   );
   
   /* 
    * Write an RGBA (RGB with alpha) color image to a png file.
    * The red, green, blue, and alpha channels should have values in [0,1].
    */
   static void write(
      const string<>&,     /* filename */
      const matrix<>&,     /* red channel */
      const matrix<>&,     /* green channel */
      const matrix<>&,     /* blue channel */
      const matrix<>&,     /* alpha channel */
      unsigned long = 16   /* bit depth (8 or 16) */
   );
};
 
} /* namespace image */
} /* namespace formats */
} /* namespace io */

#endif
