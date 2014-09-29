// Copyright (C) 2010 Charless C. Fowlkes <fowlkes@ics.uci.edu>
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

// ------------------------------------------------------------------------ 
//  June 2014, modifications by Pablo Arbelaez <arbelaez@berkeley.edu>
// ------------------------------------------------------------------------

#include "stdio.h"
#include "stdlib.h"
#include "math.h"                                       
#include "float.h"
#include "assert.h"
#include "string.h"
#include <fstream>
#include <iostream>
#include "array.hh"
#include "smatrix.hh"
#include "affinity.hh"
#include "ic.hh"
#include "mex.h"

using namespace std;


void mexFunction(int nlhs, mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    if (nrhs != 4) mexErrMsgTxt("INPUT:  (lattice.H, lattice.V, dthresh, sigma) ");
    if (nlhs != 3) mexErrMsgTxt("OUTPUT: [i,j,s] ");

    
    int dthresh;
    float sigma;
    
    dthresh = int(mxGetScalar(prhs[2]));
    sigma = float(mxGetScalar(prhs[3]));

    // copy edge info into lattice struct
    Group::DualLattice boundaries; 
    double* H = mxGetPr(prhs[0]);
    int H_h = mxGetN(prhs[0]);
    int H_w = mxGetM(prhs[0]);
    boundaries.H.resize(H_h,H_w);
    for (int i = 0; i < H_h; i++)
    {
      for (int j = 0; j < H_w; j++)
      {
        boundaries.H(i,j) = H[i*H_w+j];
      } 
    }
    double* V = mxGetPr(prhs[1]);
    int V_h = mxGetN(prhs[1]);
    int V_w = mxGetM(prhs[1]);
    boundaries.V.resize(V_h,V_w);
    for (int i = 0; i < V_h; i++)
    {
      for (int j = 0; j < V_w; j++)
      {
        boundaries.V(i,j) = V[i*V_w+j];
      } 
    }
    boundaries.width = boundaries.H.size(0);
    boundaries.height = boundaries.V.size(1);


    Group::SupportMap ic;
    Group::computeSupport(boundaries,dthresh,1.0f,ic);

    SMatrix* W = NULL;
    Group::computeAffinities2(ic,sigma,dthresh,&W);
    if (W == NULL)
    {
      mexErrMsgTxt("compute affinities failed!");
    }

  /*
   * alternate code to read in a sparse matrix from a file
   * 
    SMatrix* W;
    char *filename = static_cast<char*>(mxCalloc(mxGetN(prhs[0])+1, sizeof(char))); //mxCalloc is similar to malloc in C
    mxGetString(prhs[0],filename,mxGetN(prhs[0])+1);
    FILE* fp = fopen(filename,"r");
    if (fp != NULL)
    {
      W = new SMatrix(fp);
    }
    else
    {
      mexErrMsgTxt("Unable to open file for input.");
    }
   */

    // now pack sparse matrix into LHS.
    // compute total number of nonzero entries
    int nnz = 0;
    for (int i = 0; i < W->n; i++)
      nnz += W->nz[i];
    plhs[0] = mxCreateDoubleMatrix(nnz,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nnz,1,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(nnz,1,mxREAL);
    double* x = mxGetPr(plhs[0]); //col index
    double* y = mxGetPr(plhs[1]); //row index
    double* z = mxGetPr(plhs[2]); //values
    int ct = 0;
    for (int row = 0; row < W->n; row++) 
    {
      for (int i = 0; i < W->nz[row]; i++)
      {
        y[ct+i] = static_cast<double>(row + 1); //add one for matlab indexing
        x[ct+i] = static_cast<double>(W->col[row][i] + 1); 
        z[ct+i] = static_cast<double>(W->values[row][i]);
      } 
      ct = ct + W->nz[row];
    }
    
    delete W;
}









