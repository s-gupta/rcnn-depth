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
    if(nrhs<2)
        mexErrMsgTxt("There should be at least 2 input parameters");
    
    /* Input parameters */
    ConstMatlabMultiArray<double> pairs(prhs[0]);
    ConstMatlabMultiArray<double> base_features(prhs[1]);
    unsigned int mode = 1;
    if (nrhs>2) mode = mxGetScalar(prhs[2]);

    std::size_t n_pairs   = pairs.shape()[0];
    unsigned int n_neighs = pairs.shape()[1];
    std::size_t n_regs    = base_features.shape()[0];
    std::size_t n_feats   = base_features.shape()[1];
       
    /* Assign output parameters */
    plhs[0] = mxCreateDoubleMatrix(n_pairs,1,mxREAL);
    MatlabMultiArray<double> features(plhs[0]);
    
    if (mode==1) // end th
    {
        for (std::size_t ii=0; ii<n_pairs; ++ii)
        {
            // Get minimum th
            features[ii][0] = std::numeric_limits<double>::max();
            for (std::size_t jj=0; jj<n_neighs; ++jj)
            {
                if (pairs[ii][jj]<0)
                    break;
                features[ii][0] = std::min(features[ii][0], base_features[pairs[ii][jj]][0]);
            }
        }
    }
    else if (mode==2) // area balance
    {
        plhs[1] = mxCreateDoubleMatrix(n_pairs,1,mxREAL);
        MatlabMultiArray<double> areas(plhs[1]);
        for (std::size_t ii=0; ii<n_pairs; ++ii)
        {
            // Get area balance
            double min_area = std::numeric_limits<double>::max();
            double max_area = std::numeric_limits<double>::min();
            areas[ii][0] = 0;
            for (std::size_t jj=0; jj<n_neighs; ++jj)
            {
                if (pairs[ii][jj]<0)
                    break;
                min_area = std::min(min_area, base_features[pairs[ii][jj]][0]);
                max_area = std::max(max_area, base_features[pairs[ii][jj]][0]);
                areas[ii][0] += base_features[pairs[ii][jj]][0];
            }
            features[ii][0] = (double)min_area/(double)max_area;
        }
    }
    else if (mode==3) // Perimeters
    {
        if (nrhs<4)
            mexErrMsgTxt("For perimeter we need 4 parameters");
        ConstMatlabMultiArray<double> common_perims(prhs[3]);
        // base_features-> Perimeters

        plhs[1] = mxCreateDoubleMatrix(n_pairs,1,mxREAL);
        MatlabMultiArray<double> in_contour(plhs[1]);
        
        for (std::size_t ii=0; ii<n_pairs; ++ii)
        {
            in_contour[ii][0] = 0;
            double perim_sum = 0;
            for (std::size_t jj=0; jj<n_neighs; ++jj)
            {
                if (pairs[ii][jj]<0)
                    break;
                for (std::size_t kk=jj+1; kk<n_neighs; ++kk)
                {
                    if (pairs[ii][kk]<0)
                        break;
                    in_contour[ii][0] += common_perims[pairs[ii][jj]][pairs[ii][kk]];
                }
                perim_sum +=  base_features[pairs[ii][jj]][0];
            }
            features[ii][0] = perim_sum - 2*in_contour[ii][0];
        }
    }
    else if (mode==4) // Bounding boxes
    {
        plhs[1] = mxCreateDoubleMatrix(n_pairs,4,mxREAL);
        MatlabMultiArray<double> bboxes(plhs[1]);
        for (std::size_t ii=0; ii<n_pairs; ++ii)
        {
            bboxes[ii][0] = std::numeric_limits<double>::max();
            bboxes[ii][1] = std::numeric_limits<double>::max();
            bboxes[ii][2] = -std::numeric_limits<double>::max();
            bboxes[ii][3] = -std::numeric_limits<double>::max();
            for (std::size_t jj=0; jj<n_neighs; ++jj)
            {
                if (pairs[ii][jj]<0)
                    break;
                bboxes[ii][0] = std::min(bboxes[ii][0], base_features[pairs[ii][jj]][0]);
                bboxes[ii][1] = std::min(bboxes[ii][1], base_features[pairs[ii][jj]][1]);
                bboxes[ii][2] = std::max(bboxes[ii][2], base_features[pairs[ii][jj]][2]);
                bboxes[ii][3] = std::max(bboxes[ii][3], base_features[pairs[ii][jj]][3]);
            }
        }
    }

}
