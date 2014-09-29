// ------------------------------------------------------------------------ 
//  Copyright (C)
//  Universitat Politecnica de Catalunya BarcelonaTech (UPC) - Spain
//  University of California Berkeley (UCB) - USA
// 
//  Jordi Pont-Tuset <jordi.pont@upc.edu>
//  Pablo Arbelaez <arbelaez@berkeley.edu>
//  June 2014
// ------------------------------------------------------------------------ 
// This file is part of the MCG package presented in:
//    Arbelaez P, Pont-Tuset J, Barron J, Marques F, Malik J,
//    "Multiscale Combinatorial Grouping,"
//    Computer Vision and Pattern Recognition (CVPR) 2014.
// Please consider citing the paper if you use this code.
// ------------------------------------------------------------------------
#include "mex.h"

#include <iostream>
#include <set>
#include <map>
#include <list>
#include "matlab_multiarray.hpp"

using namespace std;

void mexFunction( int nlhs, mxArray *plhs[], 
        		  int nrhs, const mxArray*prhs[] )
{    
    /* Check for proper number of arguments */
    if (nrhs != 3) { 
    	mexErrMsgTxt("Three input arguments required."); 
    } else if (nlhs > 3) {
        mexErrMsgTxt("Too many output arguments."); 
    } 
    
    /* Input as a Multiarray */
    ConstMatlabMultiArray<double> lp(prhs[0]); /* Leaves partition */
    ConstMatlabMultiArray<double> ms(prhs[1]); /* Merging sequence */
    ConstMatlabMultiArray<double> cands(prhs[2]); /* Candidates */
    
    /* Input sizes and checks */
    size_t sm = lp.shape()[0];
    size_t sn = lp.shape()[1];
    size_t n_sons_max= ms.shape()[1]-1;
    size_t n_merges  = ms.shape()[0];
    size_t n_leaves  = ms[0][n_sons_max]-1; // n_merges+1;  --> Not valid for non-binary trees
    size_t n_regs    = n_leaves+n_merges;
    size_t n_cands   = cands.shape()[0];
    size_t n_regs_max= cands.shape()[1];
   
    
    /* Create the list of activated pixels in each leave region */
    vector<list<size_t> > x_coords(n_regs);
    vector<list<size_t> > y_coords(n_regs);
    for (size_t xx=0; xx<sm; ++xx)
    {
        for (size_t yy=0; yy<sn; ++yy)
        {
            size_t curr_leave = lp[xx][yy]-1;
            x_coords[curr_leave].push_back(xx);
            y_coords[curr_leave].push_back(yy);
        }
    }

    /* Bottom-up expansion of the lists of pixels */
    for (size_t ii=0; ii<n_merges; ++ii)
    {
        for (size_t jj=0; jj<n_sons_max; ++jj)
        {
            if (ms[ii][jj]==0)
                break;
            x_coords[n_leaves+ii].insert(x_coords[n_leaves+ii].end(),
                                         x_coords[ms[ii][jj]-1].begin(),
                                         x_coords[ms[ii][jj]-1].end());
            y_coords[n_leaves+ii].insert(y_coords[n_leaves+ii].end(),
                                         y_coords[ms[ii][jj]-1].begin(),
                                         y_coords[ms[ii][jj]-1].end());
        }
    }
    
    /* Allocate out masks */
    size_t ndim = 3;
    mwSize dims[3];
    dims[0] = sm;  dims[1] = sn;   dims[2] = n_cands;
    plhs[0] = mxCreateLogicalArray(3, dims);
    MatlabMultiArray3<bool> out_masks(plhs[0]);
    
    for (size_t ii=0; ii<n_cands; ++ii)
    {
        for (size_t jj=0; jj<n_regs_max; ++jj)
        {
            if (cands[ii][jj]==0)
                break;
            
            list<size_t>::iterator itx = x_coords[cands[ii][jj]-1].begin();
            list<size_t>::iterator ity = y_coords[cands[ii][jj]-1].begin();
            for( ; itx!=x_coords[cands[ii][jj]-1].end(); ++itx, ++ity)
            {
                out_masks[*itx][*ity][ii] = true;
            }
        }
    }
}

