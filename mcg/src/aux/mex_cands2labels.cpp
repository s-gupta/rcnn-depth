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
    if (nrhs != 2) { 
    	mexErrMsgTxt("Two input arguments required."); 
    } else if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments."); 
    } 
    
    /* Input as a Multiarray */
    ConstMatlabMultiArray<double> ms(prhs[0]); /* Merging sequence */
    ConstMatlabMultiArray<double> cands(prhs[1]); /* Candidates */
    
    /* Input sizes and checks */
    size_t n_sons_max= ms.shape()[1]-1;
    size_t n_merges  = ms.shape()[0];
    size_t n_leaves  = ms[0][n_sons_max]-1; // n_merges+1;  --> Not valid for non-binary trees
    size_t n_regs    = n_leaves+n_merges;
    size_t n_cands   = cands.shape()[0];
    size_t n_regs_max= cands.shape()[1];
   
    /* Fill leaf labels */
    vector<list<size_t> > tree_labels(n_regs);
    for (size_t ii=0; ii<n_leaves; ++ii)
        tree_labels[ii].push_back(ii);
    
    /* Bottom-up expansion of the labels */
    for (size_t ii=0; ii<n_merges; ++ii)
    {
        for (size_t jj=0; jj<n_sons_max; ++jj)
        {
            if (ms[ii][jj]==0)
                break;
            tree_labels[ms[ii][n_sons_max]-1].insert(tree_labels[ms[ii][n_sons_max]-1].end(),
                                                     tree_labels[ms[ii][jj]-1].begin(),
                                                     tree_labels[ms[ii][jj]-1].end());
        }
    }
    
    /* Create the labels for all candidates */
    vector<vector<size_t> > labels(n_cands);
    for (size_t ii=0; ii<n_cands; ++ii)
    {
        for (size_t jj=0; jj<n_regs_max; ++jj)
        {
            if (cands[ii][jj]==0)
                break;
            labels[ii].insert(labels[ii].end(),
                              tree_labels[cands[ii][jj]-1].begin(),
                              tree_labels[cands[ii][jj]-1].end());
        }
    }
    
    /* Sort all labels */
    for (size_t ii=0; ii<n_cands; ++ii)
    {
        std::sort(labels[ii].begin(),labels[ii].end());
    }
    
    /* Store labels as a cell */
    int dims_out[1]; 
    dims_out[0] = n_cands;
    plhs[0] = mxCreateCellArray(1, dims_out);
    
    for (size_t ii=0; ii<n_cands; ++ii)
    {
        mxArray* curr_entry = mxCreateNumericMatrix(1, labels[ii].size(), mxUINT32_CLASS, mxREAL);
        int* curr_pr = (int*)mxGetPr(curr_entry);
        std::size_t curr_jj = 0;
        for(vector<size_t>::const_iterator list_it=labels[ii].begin(); list_it!=labels[ii].end(); ++list_it)
        {
            curr_pr[curr_jj] = *list_it+1;
            curr_jj++;
        }
        mxSetCell(plhs[0], ii, curr_entry);
    }
}

