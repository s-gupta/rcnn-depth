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
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>

#include "smatrix.hh"
#include "ic.hh"
#include "affinity.hh"

namespace Group
{

  //////////////////////////////////////////////////////////////////////////////
  //
  // probability-of-same-segment logistic fits for color and grayscale
  // 
  // optimized for dthresh = 20
  //
  //

  float C_IC_SS(float ic)
  {
    float val = (float)1.8682;

    val += (ic / (float)0.3130)*-(float)1.3113;

    const float post = (float)1.0 / ((float)1.0 + exp(-val));
#ifdef _MSC_VER
    if (!_finite(post)) { return 0; }
#else
    if (!isfinite(post)) { return 0; }
#endif
    if (post < 0) { return 0; }
    if (post > 1) { return 1; }
    return post;
    
  }
  
  //
  // compute similarities for the set of "true" pixels in region.  
  // affinity matrix is ordered in scanline order 
  //
  void computeAffinities2(const SupportMap& icmap, const float sigma, const float dthresh, SMatrix** affinities)
  {
    int width = icmap.size(0);
    int height = icmap.size(1);

    //build a scanline order index 
    int numPixels = 0;
    Util::Array2D<int> index(width,height);
    index.init(-1);
    for (int y = 0; y < height; y++)
    {
      for (int x = 0; x < width; x++)
      {
        index(x,y) = numPixels;
        numPixels++;
      }
    }

    //sparse matrix data
    int* nz = new int[numPixels];           //number of non-zero entries in each row
    double** vals = new double*[numPixels];   //the values in each row
    int** col = new int*[numPixels];        //the column number for each value in the row  
   
    int dthreshi = (int)ceil(dthresh);
    int wd = (2*dthreshi+1);                 //window diameter 
    Util::Array1D<PointIC> connections(wd*wd);
    
    int row = 0;
    for (int x = 0; x < width; x++)
    {
        for (int y = 0; y < height; y++)
        {
            row = index(x,y);  //the row we are working on
            nz[row] = 0;       //connection count for row i
            int icIndex = 0;        //index into sparse supportMap
            for (int u = -dthreshi; u <= dthreshi; u++)
            {
              int yy = y + u;
              for (int v = -dthreshi; v <= dthreshi; v++)
              {
                int xx = x + v;
                if (xx < 0 || xx >= width) {continue;}
                if (yy < 0 || yy >= height) {continue;}
                if (u*u+v*v > dthresh*dthresh) {continue;}
              
                //increment our index into the support map
                while( icIndex < icmap(x,y).size() && 
                       icmap(x,y)(icIndex).y < yy) 
                {
                  icIndex++;
                }
                while( icIndex < icmap(x,y).size() && 
                        icmap(x,y)(icIndex).x < xx) 
                {
                  icIndex++;
                }


                float pss = 0.0;     //connection strength
                if ((u == 0) && (v == 0))
                {
                  pss = 1.0f;
                } 
                else
                {
                  float icsim = 0.0f;
                  if (icIndex < icmap(x,y).size() &&
                       icmap(x,y)(icIndex).x == xx &&
                        icmap(x,y)(icIndex).y == yy)
                  {
                    icsim = icmap(x,y)(icIndex).sim;
                    icIndex++;
                  }
                  pss = C_IC_SS(1-icsim);
                }//if (u==0) & (v==0)

                connections(nz[row]).x = xx;
                connections(nz[row]).y = yy;
                connections(nz[row]).sim = pss;
                nz[row]++;
              }//for v
            }//for u

            //fill in entries of sparse matrix
            vals[row] = new double[nz[row]];
            col[row] = new int[nz[row]];
            for (int j = 0; j < nz[row]; j++)
            {
              float val = exp( -(1-connections(j).sim) / sigma);
              assert((val >= 0.0) && (val <= 1.0));
              vals[row][j] = val;
              col[row][j] = index(connections(j).x,connections(j).y);
            }
        }//for y
    }//for x

    *affinities = new SMatrix(numPixels,nz,col,vals);
    (*affinities)->symmetrize();
  }

} //namespace Group



