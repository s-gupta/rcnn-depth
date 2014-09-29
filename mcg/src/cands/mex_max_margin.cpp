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

void mexFunction( int nlhs, mxArray *plhs[], 
        		  int nrhs, const mxArray*prhs[] )
{
    if(nrhs<3)
        mexErrMsgTxt("There should be at least 3 input parameters");
    
    /* Input parameters */
    ConstMatlabMultiArray<double> pairs(prhs[0]);
    ConstMatlabMultiArray<double> scores(prhs[1]);
    ConstMatlabMultiArray<double> base_intersections(prhs[2]);
    double theta = 0.75;
    if (nrhs>3) theta = mxGetScalar(prhs[3]);

    std::size_t n_pairs   = pairs.shape()[0];
    unsigned int n_neighs = pairs.shape()[1];
    std::size_t n_regs    = base_intersections.shape()[0];
    
    /* Output allocation */
    plhs[0] = mxCreateDoubleMatrix(n_pairs,1,mxREAL);
    MatlabMultiArray<double> sel_ids(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(n_pairs,1,mxREAL);
    MatlabMultiArray<double> new_scores(plhs[1]);
    
    /* Current maximum overlap */
    std::vector<double> max_J(n_pairs,0);
    
    /* Remaining candidates */
    std::set<std::size_t> remaining;
    for (std::size_t ii=0; ii<n_pairs; ++ii)
        remaining.insert(ii);
    
    /**/
    std::size_t n_cands(0);
    for (std::size_t ii=0; ii<n_pairs; ++ii)
    {
        /* Recompute remaining scores */
        std::set<std::size_t>::iterator it = remaining.begin();
        double curr_max = -std::numeric_limits<double>::max();
        std::size_t next_to_add;
//         std::cout << "Remain: ";
        for( ; it!=remaining.end(); ++it)
        {
//             std::cout << *it+1 << ", ";
            double curr_score = theta*scores[*it][0] - (1-theta)*max_J[*it];
            if (curr_score>curr_max)
            {
                curr_max = curr_score;
                next_to_add = *it;
            }
        }
//         std::cout << std::endl;
//         std::cout << next_to_add+1 << std::endl;
        
        /* Add the candidate */
        remaining.erase(next_to_add);
        sel_ids[n_cands][0] = next_to_add + 1;
        new_scores[n_cands][0] = curr_max;
        n_cands++;
        
        /* Update max overlap with the new candidate added */
        for(it = remaining.begin(); it!=remaining.end(); ++it)
        {
            /* Compute overlap between ramaining and next_to_add */
            double curr_inters = 0;
            double area1 = 0;
            double area2 = 0;
            for (std::size_t jj=0; jj<n_neighs; ++jj)
            {
                if (pairs[next_to_add][jj]<0)
                    break;
                for (std::size_t kk=0; kk<n_neighs; ++kk)
                {
                    if (pairs[*it][kk]<0)
                        break;
                    curr_inters += base_intersections[pairs[next_to_add][jj]][pairs[*it][kk]];
                }
                area1 += base_intersections[pairs[next_to_add][jj]][pairs[next_to_add][jj]];
            }
            for (std::size_t kk=0; kk<n_neighs; ++kk)
            {
                if (pairs[*it][kk]<0)
                    break;
                area2 += base_intersections[pairs[*it][kk]][pairs[*it][kk]];
            }

            double curr_J = (double)curr_inters/(double)(area1+area2-curr_inters);
 
//             std::cout << next_to_add+1 << " / " << *it+1 << ": " << area1 << ", " << area2 << ", " << curr_inters << " || " << curr_J  << std::endl;
                    
            /* Is it higher than the previous one? */
            max_J[*it] = std::max(max_J[*it], curr_J);
        }
    }
}
