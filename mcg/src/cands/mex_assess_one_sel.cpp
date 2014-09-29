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
// ------------------------------------------------------------------------#include "mex.h"
#include "matlab_multiarray.hpp"
#include <iostream>
#include <list>
#include <set>
#include <algorithm>

void mexFunction( int nlhs, mxArray *plhs[], 
        		  int nrhs, const mxArray*prhs[] )
{   
    if(nrhs!=3)
        mexErrMsgTxt("There should be 3 input parameters");
    
    /* Input parameters */
    ConstMatlabMultiArray<double> tps(prhs[0]);
    ConstMatlabMultiArray<double> fps(prhs[1]);
    ConstMatlabMultiArray<double> pairs(prhs[2]);
    
    std::size_t n_pairs   = pairs.shape()[0];
    std::size_t n_neighs  = pairs.shape()[1];
       
    /* Assign output parameters */
    plhs[0] = mxCreateDoubleMatrix(n_pairs,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(n_pairs,1,mxREAL);

    MatlabMultiArray<double> pairs_tps(plhs[0]);
    MatlabMultiArray<double> pairs_fps(plhs[1]);
    
    for(std::size_t ii=0; ii<n_pairs; ++ii)
    {
        pairs_tps[ii][0] = 0;
        pairs_fps[ii][0] = 0;
        for(std::size_t jj=0; jj<n_neighs; ++jj)
        {
            if (pairs[ii][jj]==0)
                break;
            pairs_tps[ii][0] += tps[0][pairs[ii][jj]-1];
            pairs_fps[ii][0] += fps[0][pairs[ii][jj]-1];
        }
    }
}

