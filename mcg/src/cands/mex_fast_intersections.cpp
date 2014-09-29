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
#include "matlab_multiarray.hpp"
#include <iostream>
#include <list>
#include <set>
#include <algorithm>

void vector_union(std::vector<unsigned int>& vec1, std::vector<unsigned int>& vec2, std::vector<unsigned int>& vec_union)
{
    vec_union.resize(vec1.size()+vec2.size());
    std::vector<unsigned int>::iterator it = std::set_union(
            vec1.begin(), vec1.end(),
            vec2.begin(), vec2.end(), vec_union.begin());
    vec_union.resize(it-vec_union.begin());
}

void vector_intersection(std::vector<unsigned int>& vec1, std::vector<unsigned int>& vec2, std::vector<unsigned int>& vec_inters)
{
    vec_inters.resize(vec1.size()+vec2.size());
    std::vector<unsigned int>::iterator it = std::set_intersection(
            vec1.begin(), vec1.end(),
            vec2.begin(), vec2.end(), vec_inters.begin());
    vec_inters.resize(it-vec_inters.begin());
}

void mexFunction( int nlhs, mxArray *plhs[], 
        		  int nrhs, const mxArray*prhs[] )
{   
    if(nrhs<3)
        mexErrMsgTxt("There should be at least 3 input parameters");
    
    /* Input parameters */
    ConstMatlabMultiArray<double> leaves_part(prhs[0]);
    ConstMatlabMultiArray<double> mer_seq    (prhs[1]);
    ConstMatlabMultiArray<double> base_areas (prhs[2]);

    std::size_t n_regs    = base_areas.shape()[0];

    /* Get leaves of original regions */
    std::size_t n_merges = mer_seq.shape()[0];
    std::size_t n_sons_max= mer_seq.shape()[1]-1;
    std::size_t n_leaves = n_regs-n_merges;
    std::vector<std::vector<unsigned int> > leaves(n_regs);
    for (std::size_t ii=0; ii<n_leaves; ++ii)
        leaves[ii].push_back(ii);
    for (std::size_t ii=0; ii<n_merges; ++ii)
    {
        unsigned int parent = mer_seq[ii][n_sons_max];
        for (std::size_t jj=0; jj<n_sons_max; ++jj)
        {
            if (mer_seq[ii][jj]<0)
                break;
        
            unsigned int son = mer_seq[ii][jj];
            leaves[parent].insert(leaves[parent].begin(),leaves[son].begin(),leaves[son].end());
        }
        std::sort(leaves[parent].begin(), leaves[parent].end());
    }
    
    /* Assign output parameters */
    plhs[0] = mxCreateDoubleMatrix(n_regs,n_regs,mxREAL);
    MatlabMultiArray<double> intersections(plhs[0]);
    
    for (std::size_t ii=0; ii<n_regs; ++ii)
    {
        for (std::size_t jj=ii+1; jj<n_regs; ++jj)
        {
            intersections[ii][jj] = 0;
            intersections[jj][ii] = 0;
            if (jj>=n_leaves) // Leaves do no intersect between them
            {
                std::vector<unsigned int> curr_inters;
                vector_intersection(leaves[ii],leaves[jj],curr_inters);
                for  (std::size_t kk=0;kk<curr_inters.size(); ++kk)
                {
                    intersections[ii][jj] += base_areas[curr_inters[kk]][0];
                    intersections[jj][ii] += base_areas[curr_inters[kk]][0];
                }
            }
        }
    }    
}
