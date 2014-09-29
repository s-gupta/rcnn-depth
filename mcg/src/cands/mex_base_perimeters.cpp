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
#include <map>
#include <algorithm>

using namespace std;

void mexFunction( int nlhs, mxArray *plhs[], 
        		  int nrhs, const mxArray*prhs[] )
{
    if(nrhs!=6)
        mexErrMsgTxt("There should be 6 input parameters");
    
    /* Input parameters */
    double nregs = mxGetScalar(prhs[0]);
    ConstMatlabMultiArray<double> pairs_min(prhs[1]);
    ConstMatlabMultiArray<double> pairs_max(prhs[2]);
    ConstMatlabMultiArray<double> pairs_ids(prhs[3]);
    ConstMatlabMultiArray<double> pairs_matrix(prhs[4]);
    ConstMatlabMultiArray<double> feats(prhs[5]);
    
    size_t sx   = pairs_matrix.shape()[0];
    size_t sy   = pairs_matrix.shape()[1];
    size_t n_neighs = max(pairs_ids.shape()[0],pairs_ids.shape()[1]);
    
    /*-------------------------------------------------------------*/
    /*  */
    /*-------------------------------------------------------------*/
    map<double,double> lut;
    for (size_t ii=0; ii<n_neighs; ++ii)
        lut.insert(pair<double,double>(pairs_ids[ii][0],ii));
    
    vector<double> counts(n_neighs,0);
        
    /* Scan horizontal contours*/
    for (std::size_t xx=0; xx<sx; xx+=2)
        for (std::size_t yy=1; yy<sy; yy+=2)
            if (pairs_matrix[xx][yy]>0)
            {
                map<double,double>::iterator it = lut.find(pairs_matrix[xx][yy]);
                counts[it->second] += feats[xx][yy];
            }

    /* Scan vertical contours*/
    for (std::size_t xx=1; xx<sx; xx+=2)
        for (std::size_t yy=0; yy<sy; yy+=2)
            if (pairs_matrix[xx][yy]>0)
            {
                map<double,double>::iterator it = lut.find(pairs_matrix[xx][yy]);
                counts[it->second] += feats[xx][yy];
            }
    /*-------------------------------------------------------------*/
 
    
    /* Output allocation */
    plhs[0] = mxCreateDoubleMatrix(nregs,nregs,mxREAL);
    MatlabMultiArray<double> perims(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(nregs,1,mxREAL);
    MatlabMultiArray<double> border_perims(plhs[1]);
    
    for (size_t ii=0; ii<n_neighs; ++ii)
    {
        if (pairs_min[ii][0]==0)
            border_perims[pairs_max[ii][0]-1][0] = counts[ii];
        else
            perims[pairs_min[ii][0]-1][pairs_max[ii][0]-1] = counts[ii];
    }
}
    
