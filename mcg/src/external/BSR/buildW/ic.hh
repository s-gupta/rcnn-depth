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


#ifndef IC_HH
#define IC_HH

#include "image.hh"

//
// intervening contour computation
//
namespace Group
{
  // data types for storing sparse ic info
  struct DualLattice
  {
    Util::Image H;
    Util::Image V;
    int width;
    int height;
  };
 


  //
  // the (1-ic) to a given point
  //
  struct PointIC
  {
    int x, y;
    float sim;
  };

  //
  // (1-ic) from each pixel to some set of neighbors
  //
  typedef Util::Array2D< Util::Array1D<PointIC> > SupportMap;

  //
  // given a pb image and window radius, computes a support map (1-ic) for each
  // pixel out to the given radius.
  //
  void computeSupport(const DualLattice& boundaries, const int wr, const float thresh, SupportMap& support);
 
  //
  // given a pb image, a pixel (x0,y0), and a pb threshold, compute the intervening-contour
  // weight to all pixels inside a given box of radius "wr" subject to the threshold "thresh".
  // pb is max-accumulated along bresenham lines.  the result is stored in scanline order as
  // an array of <PointIC>.
  //
  void interveningContour(const DualLattice& boundaries, const float thresh,
                          const int x0, const int y0, const int wr,
                          Util::Array1D<PointIC> &adj, int &count);

  //
  // compute (1 - max over lattice energies on a straightline path connecting p1 and p2)
  //
  void interveningContour (const DualLattice& boundaries, const int x1, const int y1, 
                           const int x2, const int y2, float& icsim);

} //namespace Group

#endif

