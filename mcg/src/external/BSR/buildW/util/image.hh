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

#ifndef IMAGE_HH
#define IMAGE_HH

#include "array.hh"

namespace Util 
{
  typedef Array2D<float> Image;
  typedef Array3D<float> ImageStack;
  typedef Array2D<float> Matrix;


  enum RGB_CHANNELS {RGB_R=0,RGB_G=1,RGB_B=2};
  enum LAB_CHANNELS {LAB_L=0,LAB_A=1,LAB_B=2};

  //
  // read in a jpeg file into an ImageStack.
  // if the jpeg file is grayscale then the resulting ImageStack only has 1 layer
  // otherwise it has 3 layers corresponding to the RGB colorspace.
  //
  bool readJpegFile (const char *filespec, ImageStack& im);

  //
  // write out a grayscale image to a jpeg file.
  // if normalize=true, then the range of the image is adjusted to use the full scale
  // if jet=true then the image is written in pseudocolor rather than grayscale
  //
  bool writeJpegFile (const Image& im, const char *filespec, 
                      const bool normalize = true, const bool jet = false);

  //
  // convert an RGB imagestack into an 1976 CIE L*a*b* imagestack.
  //
  void rgb2lab(const ImageStack& rgb, ImageStack& lab);

  //
  // normalize a given Lab image stack so that values all lie in [0,1]
  //
  void labNormalize(ImageStack& lab);

  //
  // create a translated version of this image where
  // the old image appears embedded in a new image of
  // size [newwidth x newheight].  undefined pixels 
  // are filled in with value fill.
  //
  void getTranslated(const Image& im, const int xoffset, const int yoffset,
                     const int newwidth, const int newheight, 
                     const float fill, Image& translated);

  //
  // create a resized version of this image of size [newwidth x newheight]
  // if bilinear = true, use bilinear interpolation
  // otherwise use bicubic interpolation
  //
  void getScaled (const Image& im, const int newwidth, const int newheight, 
                  const bool bilinear, Image& scaled);

  //
  // create a rotated version of this image
  // using an appropriate affine transform.
  // if bilinear = true, use bilinear interpolation
  // otherwise use bicubic.
  //
  void getRotated (const Image& im, const float theta, 
                   const bool bilinear, Image& rotated);


  //
  // returns a new image which is an affine
  // transformed version of this image.
  // newimage = A*image.  the new image
  // is of size (height, width) such that
  // the corners of the old image are transformed
  // to locations inside the new image
  //
  // if bilinear is TRUE then use bilinear interpolation
  // otherwise use bicubic B-spline interpolation
  //
  void getTransformed (const Image& im, const Matrix& A, 
                       const int width, const int height, 
                       const int xoffset, const int yoffset, 
                       const bool bilinear, Image& transformed);

  //
  // filters the image via convolution with the given
  // kernel and returns the resulting image.  kernel
  // must have odd dimensions.
  //
  void getFiltered (const Image& im, const Image& kernel, Image& filtered);

  //
  // filter at the given radius.
  // the resulting value at a pixel is the n'th-order value
  // in a circular window cetnered at the pixel
  // 0 <= order <= 1, 0 is smallest value, 1 is largest.
  // Call with 0.5 to get median filtering.
  //
  void getPctFiltered (const Image& im, const float radius, 
                       const float order, Image& filtered);


  //
  // filter at the given radius.
  // the resulting value at a pixel is the maximum value
  // in a circular window cetnered at the pixel
  // 
  void getMaxFiltered (const Image& im, const float radius,
                       Image& filtered);

}
#endif

