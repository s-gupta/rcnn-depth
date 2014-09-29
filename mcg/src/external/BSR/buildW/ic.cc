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


#include <assert.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <ctype.h>
#include <memory.h>

#include "image.hh"
#include "array.hh"

#include "ic.hh"

namespace Group
{
  //
  // given a pb image and window radius, computes a support map for each
  // pixel out to the given radius.
  //
  void computeSupport(const DualLattice& boundaries, const int wr,
                       const float thresh, SupportMap& support)
  {
    support.resize(boundaries.width,boundaries.height);
    Util::Array1D<PointIC> adj;
    int count = 0;

    //Util::Message::startBlock(boundaries.width,"computing support");
    for (int x = 0; x < boundaries.width; x++)
    {
      //Util::Message::stepBlock();
      for (int y = 0; y < boundaries.height; y++)
      {
        interveningContour(boundaries,thresh,x,y,wr,adj,count);    
        Util::Array1D<PointIC> map(count);
        for (int i = 0; i < count; i++)
        {
          map(i) = adj(i);    
          const int ix = map(i).x;
          const int iy = map(i).y;
          assert(ix >= 0);
          assert(ix < boundaries.width);
          assert(iy >= 0);
          assert(iy < boundaries.height);
        }
        support(x,y) = map;
      }
    }
    //Util::Message::endBlock();
  }

  //
  // Walk the bresenham line from (x0,y0) to (x2,y2) ignoring any points outside
  // a circular window of radius wr. 
  // (x1,y1) and (x3,y3) should be on either side of (x2,y2)
  //
  // For each line, stop if the max-accumulated pb is > thresh. 
  //
  // Append any points on the line (for which the line is the best approximant
  // and distance from x0 is less than wr) to scanlines array.
  //
  // points is preallocated scratch space.
  // scanCount and scanLines store the results
  //
  void
  ic_walk (const DualLattice& boundaries, const float thresh, 
           const int x0, const int y0, 
           const int x1, const int y1,
           const int x2, const int y2,
           const int x3, const int y3,
           const int wr,
           Util::Array1D <PointIC> &points, 
           Util::Array1D <int>&scanCount,
           Util::Array2D <PointIC> &scanLines)
  {
      const int width = boundaries.width;
      const int height = boundaries.height;

      // the predicate that uses long longs will overflow if the image
      // is too large
      assert (2*wr+1 < 1000);

      // make sure points array is big enough
      assert ((int) points.size () >= 4*wr+2);

      // make sure scan arrays are the right size
      assert ((int)scanCount.size() == 2*wr+1);
      assert ((int)scanLines.size(0) == 2*wr+1);
      assert ((int)scanLines.size(1) == 2*wr+1);

      //sanity check the points
      assert (x0 >= 0 && x0 < width);
      assert (y0 >= 0 && y0 < height);
      assert (x1 >= 0 && x1 < width);
      assert (y1 >= 0 && y1 < height);
      assert (x2 >= 0 && x2 < width);
      assert (y2 >= 0 && y2 < height);
      assert (x3 >= 0 && x3 < width);
      assert (y3 >= 0 && y3 < height);

      // make sure points are all distinct
      assert (x0 != x1 || y0 != y1);
      assert (x0 != x2 || y0 != y2);
      assert (x0 != x3 || y0 != y3);
      assert (x1 != x2 || y1 != y2);
      assert (x1 != x3 || y1 != y3);
      assert (x2 != x3 || y2 != y3);

      //constants used in testing whether this is
      //the best path
      const long long dx1 = x1 - x0; 
      const long long dy1 = y1 - y0;
      const long long dx2 = x2 - x0;
      const long long dy2 = y2 - y0;
      const long long dx3 = x3 - x0;
      const long long dy3 = y3 - y0;
      const long long dot11 = dx1 * dx1 + dy1 * dy1;
      const long long dot22 = dx2 * dx2 + dy2 * dy2;
      const long long dot33 = dx3 * dx3 + dy3 * dy3;

      // compute dx,dy for the bresenham line
      const int dx = x2 - x0;
      const int dy = y2 - y0;
      const int adx = abs (dx);
      const int ady = abs (dy);

      // figure out what octant we're in for the bresenham algorithm;
      // octant i covers pi/4 * [i,i+1)
      int octant = -1;
      if (dx > 0 && dy >= 0) {           // quadrant 0
        octant = (adx > ady) ? 0 : 1;
      } else if (dx <= 0 && dy > 0) {    // quadrant 1
        octant = (adx < ady) ? 2 : 3;
      } else if (dy <= 0 && dx < 0) {    // quadrant 2
        octant = (adx > ady) ? 4 : 5;
      } else if (dx >= 0 && dy < 0) {    // quadrant 3
        octant = (adx < ady) ? 6 : 7;
      } else {
        assert (0);
      }

      // t is our bresenham counter
      int t = 0;
      switch (octant)
      {
        case 0: t = -adx; break;
        case 1: t = -ady; break;
        case 2: t = -ady; break;
        case 3: t = -adx; break;
        case 4: t = -adx; break;
        case 5: t = -ady; break;
        case 6: t = -ady; break;
        case 7: t = -adx; break;
        default: assert (0);
      }

      // maxpb contains the max-accumulation of pb from (x0,y0) to (x,y)
      // on the bresenham line. 
      float maxpb = 0.0f; 

      // (xi,yi) is our current location on the bresenham line
      int xi = x0;
      int yi = y0;

      // accumulate the points in the order we find them
      int count = 0;
      int oldx = xi;
      int oldy = yi;

      //walk the line
      while (xi != x2 || yi != y2)
      {
        // step one pixel on the bresenham line
        switch (octant)
        {
          case 0:
            xi++; t += (ady << 1);
            if (t > 0) { yi++; t -= (adx << 1); }
            break;
          case 1:
            yi++; t += (adx << 1);
            if (t > 0) { xi++; t -= (ady << 1); }
            break;
          case 2:
            yi++; t += (adx << 1);
            if (t > 0) { xi--; t -= (ady << 1); }
            break;
          case 3:
            xi--; t += (ady << 1);
            if (t > 0) { yi++; t -= (adx << 1); }
            break;
          case 4:
            xi--; t += (ady << 1);
            if (t > 0) { yi--; t -= (adx << 1); }
            break;
          case 5:
            yi--; t += (adx << 1);
            if (t > 0) { xi--; t -= (ady << 1); }
            break;
          case 6:
            yi--; t += (adx << 1);
            if (t > 0) { xi++; t -= (ady << 1); }
            break;
          case 7:
            xi++; t += (ady << 1);
            if (t > 0) { yi--; t -= (adx << 1); }
            break;
          default:
            assert (0);
        }

        // Figure out if the bresenham line from (x0,y0) to (x2,y2) is the
        // best approximant we will see for the line from (x0,y0) to (xi,yi).
        // We need:
        //              T(i,2) < T(i,1) && T(i,2) <= T(i,3)
        // Where T(a,b) is the angle between the two lines (x0,y0)-(xa,ya)
        // and (x0,y0)-(xb,yb).  
        // We can compute an exact integer predicate; let C be the square
        // of the cosine of T:
        //              C(i,2) > C(i,1) && C(i,2) >= C(i,3)
        // Use the identity:
        //              cos(t) = a.b/|a||b|
        // Square and cross-multiply to get rid of the divides and square
        // roots.
        // Note that we first check to see if T(i,2) == 0, in which case 
        // the line is a perfect approximant.

        const long long dxi = xi - x0;
        const long long dyi = yi - y0;
        const long long dotii = dxi * dxi + dyi * dyi;
        const long long doti1 = dxi * dx1 + dyi * dy1;
        const long long doti2 = dxi * dx2 + dyi * dy2;
        const long long doti3 = dxi * dx3 + dyi * dy3;


        const bool good = (doti2*doti2 == dotii*dot22)
                          || (dot11*doti2*doti2 > dot22*doti1*doti1
                              && dot33*doti2*doti2 >= dot22*doti3*doti3);


        // otherwise accumulate the pb value if we've crossed an edge
        float intersected = 0.0f;
        if (oldx == xi)
        {
          if (yi > oldy)
          {
            intersected = boundaries.H(xi,yi);  
          }
          else if (yi < oldy)
          {
            intersected = boundaries.H(oldx,oldy);  
          }
        } 
        else if (oldy == yi)
        {
          if (xi > oldx)
          {
            intersected = boundaries.V(xi,yi);  
          }
          else if (xi < oldx)
          {
            intersected = boundaries.V(oldx,oldy);  
          }
        }
        else
        {
          if ((xi > oldx) && (yi > oldy))  //down to right
          {
            intersected = std::max(boundaries.H(oldx,yi),intersected); 
            intersected = std::max(boundaries.H(xi,yi),intersected); 
            intersected = std::max(boundaries.V(xi,oldy),intersected); 
            intersected = std::max(boundaries.V(xi,yi),intersected); 
          }
          else if ((xi > oldx) && (yi < oldy)) //up to right
          {
            intersected = std::max(boundaries.H(oldx,oldy),intersected); 
            intersected = std::max(boundaries.H(xi,oldy),intersected); 
            intersected = std::max(boundaries.V(xi,oldy),intersected); 
            intersected = std::max(boundaries.V(xi,yi),intersected); 
          }
          else if ((xi < oldx) && (yi > oldy)) //down to left
          {
            intersected = std::max(boundaries.H(oldx,yi),intersected); 
            intersected = std::max(boundaries.H(xi,yi),intersected); 
            intersected = std::max(boundaries.V(oldx,oldy),intersected); 
            intersected = std::max(boundaries.V(oldx,yi),intersected); 
          }
          else if ((xi < oldx) && (yi < oldy)) //up to left
          {
            intersected = std::max(boundaries.H(oldx,oldy),intersected); 
            intersected = std::max(boundaries.H(xi,oldy),intersected); 
            intersected = std::max(boundaries.V(oldx,oldy),intersected); 
            intersected = std::max(boundaries.V(oldx,yi),intersected); 
          }
        }
        maxpb = std::max(maxpb,intersected);
        oldx = xi;
        oldy = yi;

        // if the approximation is not good, then skip this point
        if (!good) { continue; }

        // if the accumulated pb is too high, then stop
        if (maxpb > thresh)
        {
          break;
        }

        // record this connection
        PointIC p;
        p.x = xi;
        p.y = yi;
        p.sim = 1.0f - maxpb;
        points(count) = p;
        count++;
      }

      // add our list of points to scanLines; we have to reverse
      // the order in octants 2,3,4,5
      switch (octant)
      {
        case 0:
        case 1:
        case 6:
        case 7:
          for (int i = 0; i < count; i++)
          {
            const int yind = points(i).y - y0 + wr;
            scanLines(yind,scanCount(yind)++) = points (i);
          }
          break;
        case 2:
        case 3:
        case 4:
        case 5:
          for (int i = count - 1; i >= 0; i--)
          {
            const int yind = points(i).y - y0 + wr;
            scanLines(yind,scanCount(yind)++) = points (i);
          }
          break;
        default:
            assert (0);
      }
  }


  //
  // given a pb image, a pixel (x0,y0), and a pb threshold, compute the intervening-contour 
  // weight to all pixels inside a given box of radius "wr" subject to the threshold "thresh".
  // pb is max-accumulated along bresenham lines.  the result is stored in scanline order as
  // a list of points and their pb value.
  //
  void interveningContour(const DualLattice& boundaries, const float thresh,
                          const int x0, const int y0, const int wr,
                          Util::Array1D<PointIC> &adj, int &count)
  {
      const int width = boundaries.width;
      const int height = boundaries.height;

      // make sure (x0,y0) is valid
      assert (x0 >= 0 && x0 < width);
      assert (y0 >= 0 && y0 < height);

      // make sure adj array is big enough
      adj.resize((2*wr+1)*(2*wr+1));

      // allocate space for lists of pixels; this operation is O(1)
      // since the space need not be initialized.
      Util::Array2D <PointIC>scanLines(2*wr+1,2*wr+1);

      // we need to keep track of the length of the scan lines
      Util::Array1D <int>scanCount(2*wr+1);
      scanCount.init(0);

      // scratch space for ic_walk() function
      Util::Array1D<PointIC> scratch(4*wr+2);

      // the rectangle of interest, a square with edge of length
      // 2*wr+1 clipped to the image dimensions
      const int rxa = std::max(0,x0-wr);
      const int rya = std::max(0,y0-wr);
      const int rxb = std::min(x0+wr,width-1);
      const int ryb = std::min(y0+wr,height-1);

      // walk around the boundary, collecting points in the scanline array
      // first walk around the rectangle boundary clockwise for theta = [pi,0]
      //std::cerr << "[" << x0 << "," << y0 << "]  ";
      //std::cerr << "(" << rxa << "," << rya << ")-(" << rxb << "," << ryb << ")" << std::endl;
      if (x0 > rxa) // left 
      {
        if ((y0 > 0) && (y0 < ryb))
        {
          ic_walk(boundaries, thresh, x0,y0, rxa,y0-1, rxa,y0, rxa,y0+1, wr, scratch,scanCount,scanLines);
        }
        for (int y = y0-1; y > rya; y--)
        {
          ic_walk(boundaries, thresh, x0,y0, rxa,y-1, rxa,y, rxa,y+1, wr, scratch,scanCount,scanLines);
        }
      }

      if (x0 > rxa+1 || y0 > rya+1 || ((x0 > rxa) && (y0 > rya)) ) // top-left
      {
        ic_walk(boundaries, thresh, x0,y0, rxa,rya+1, rxa,rya, rxa+1,rya, wr, scratch, scanCount, scanLines);
      }
      if ( ((x0 == rxa) && (y0 == rya+1)) || ((x0 == rxa+1) && (y0 == rya)) )
      {
        PointIC pnt;
        pnt.x = rxa;
        pnt.y = rya;
        pnt.sim = 1.0f;
        const int yind = pnt.y - y0 + wr;
        scanLines(yind,scanCount(yind)++) = pnt;
      }

      if (y0 > rya) // top
      {
        for (int x = rxa+1; x < rxb; x++)
        {
          ic_walk(boundaries, thresh, x0,y0, x-1,rya, x,rya, x+1,rya, wr, scratch, scanCount, scanLines);
        }
      }

      if ((y0 > rya+1) || (x0 < rxb-1) || ((y0 > rya) && (x0 < rxb)) ) // top-right
      {
        ic_walk(boundaries, thresh, x0,y0, rxb-1,rya, rxb,rya, rxb,rya+1, wr, scratch, scanCount, scanLines);
      }
      if ( ((x0 == rxb-1) && (y0 == rya)) || ((x0 == rxb) && (y0 == rya+1)) )
      {
        PointIC pnt;
        pnt.x = rxb;
        pnt.y = rya;
        pnt.sim = 1.0f;
        const int yind = pnt.y - y0 + wr;
        scanLines(yind,scanCount(yind)++) = pnt;
      }


      if (x0 < rxb) // right
      {
        for (int y = rya+1; y < y0; y++)
        {
          ic_walk(boundaries, thresh, x0,y0, rxb,y-1, rxb,y, rxb,y+1, wr, scratch, scanCount, scanLines);
        }
      }

      // now counterclockwise for theta = (pi,0)
      if (x0 > rxa) // left
      {
        for (int y = y0+1; y < ryb; y++)
        {
          ic_walk(boundaries, thresh, x0,y0, rxa,y-1, rxa,y, rxa,y+1, wr, scratch, scanCount, scanLines);
        }
      }

      if ((x0 > rxa+1) || (y0 < ryb-1) || ((x0 > rxa) && (y0 < ryb))) // bottom-left
      {
        ic_walk(boundaries, thresh, x0,y0, rxa,ryb-1, rxa,ryb, rxa+1,ryb, wr, scratch, scanCount, scanLines);
      }
      if ( ((x0 == rxa) && (y0 == ryb-1)) || ((x0 == rxa+1) && (y0 == ryb)) )
      {
        PointIC pnt;
        pnt.x = rxa;
        pnt.y = ryb;
        pnt.sim = 1.0f;
        const int yind = pnt.y - y0 + wr;
        scanLines(yind,scanCount(yind)++) = pnt;
      }
      if (y0 < ryb) // bottom
      {
        for (int x = rxa+1; x < rxb; x++)
        {
          ic_walk(boundaries, thresh, x0,y0, x-1,ryb, x,ryb, x+1,ryb, wr, scratch, scanCount, scanLines);
        }
      }

      if ((y0 < ryb-1) || (x0 < rxb-1) || ((y0 < ryb) && (x0 < rxb))) // bottom-right
      {
        ic_walk(boundaries, thresh, x0,y0, rxb-1,ryb, rxb,ryb, rxb,ryb-1, wr, scratch, scanCount, scanLines);
      }
      if ( ((x0 == rxb-1) && (y0 == ryb)) || ((x0 == rxb) && (y0 == ryb-1)) )
      {
        PointIC pnt;
        pnt.x = rxb;
        pnt.y = ryb;
        pnt.sim = 1.0f;
        const int yind = pnt.y - y0 + wr;
        scanLines(yind,scanCount(yind)++) = pnt;
      }

      if (x0 < rxb) // right
      {
        for (int y = ryb-1; y > y0; y--)
        {
          ic_walk(boundaries, thresh, x0,y0, rxb,y-1, rxb,y, rxb,y+1, wr, scratch, scanCount, scanLines);
        }
        if ((y0 > 0) && (y0 < ryb))
        {
          ic_walk(boundaries, thresh, x0,y0, rxb,y0-1, rxb,y0, rxb,y0+1, wr, scratch, scanCount, scanLines);
        }
      }


      for (int y = 0; y < 2*wr+1; y++)
      {
        int len = scanCount(y);
        assert (len >= 0 && len <= 2*wr+1);

        if (y + y0 - wr < 0)
        {
          assert(len == 0);    
        }

        // check that pixels are in the right row
        for (int i = 0; i < len; i++)
        {
          assert(scanLines(y,i).y == y+y0-wr);
        }

        // check that pixels in each row are in increasing order
        for (int i = 0; i < len-1; i++)
        {
          assert(scanLines(y,i).x < scanLines(y,i+1).x);
        }
      }

      // construct the adjacency list
      count = 0;
      for (int y = 0; y < 2*wr+1; y++)
      {
        int len = scanCount(y);
        for (int i = 0; i < len; i++)
        {
          adj(count++) = scanLines(y,i);
        }
      }
  }


  // compute the max over lattice energies on a straightline 
  // path connecting p1 and p2
/*
  void interveningContour (const DualLattice& boundaries, const int x1, const int y1, 
                            const int x2, const int y2, float& icsim)
  {
      float maxpb = 0;

      const int width = boundaries.pb.size(0);
      const int height = boundaries.pb.size(1);
      const int dx = x2 - x1;
      const int dy = y2 - y1;
      const int steps = std::max (abs (dx), abs (dy));
      if (steps == 0) { return; }

      const float xincr = (float) dx / (float) steps;
      const float yincr = (float) dy / (float) steps;

      float x = x1;
      float y = y1;
      float olddist = boundaries.dist(x1,y1);
      for (int k = 0; k < steps; k++)
      {
        x += xincr;
        y += yincr;
        const int xi = (int) rint (x);
        const int yi = (int) rint (y);
        float newdist = boundaries.dist(xi,yi);
        if ( ((olddist >= 0) && (newdist < 0)) || 
              ((olddist <= 0) && (newdist > 0)) )
        {
          maxpb = std::max(maxpb, boundaries.pb(xi,yi));
        }
        olddist = newdist;
      }
      icsim = 1-maxpb;
  }
*/
} //namespace Group
