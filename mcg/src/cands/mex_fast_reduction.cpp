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
    ConstMatlabMultiArray<double> base_areas(prhs[1]);
    ConstMatlabMultiArray<double> base_intersections(prhs[2]);
    double J_th = 1;
    if (nrhs>3) J_th = mxGetScalar(prhs[3]);

    std::size_t n_pairs   = pairs.shape()[0];
    unsigned int n_neighs = pairs.shape()[1];
    std::size_t n_regs    = base_areas.shape()[0];
    
    /* Stack all candidates */
    std::list<std::vector<int> > all_pairs;
    for (std::size_t ii=0; ii<n_pairs; ++ii)
    {
        bool found = 0;
        std::list<std::vector<int> >::iterator it = all_pairs.begin();
//         std::cout << "Pair " << ii+1  << std::endl;
        for (; it!=all_pairs.end(); ++it)
        {
            double curr_inters = 0;
            double area1 = 0;
            double area2 = 0;
            for (std::size_t jj=0; jj<it->size(); ++jj)
            {
                if ((*it)[jj]<0)
                    break;
                for (std::size_t kk=0; kk<pairs[ii].size(); ++kk)
                {
                    if (pairs[ii][kk]<0)
                        break;
//                     std::cout << (*it)[jj]+1 << ", " <<  pairs[ii][kk]+1 << ": "<< base_intersections[(*it)[jj]][pairs[ii][kk]] << std::endl;
                    curr_inters += base_intersections[(*it)[jj]][pairs[ii][kk]];
                }
                area1 += base_areas[(*it)[jj]][0];
            }
            for (std::size_t kk=0; kk<pairs[ii].size(); ++kk)
            {
                if (pairs[ii][kk]<0)
                    break;
                area2 += base_areas[pairs[ii][kk]][0];
            }

            double curr_J = (double)curr_inters/(double)(area1+area2-curr_inters);
            if (curr_J>=J_th)
            {
//                 std::cout << curr_J << " || " <<  ii << "/";
//                 for (std::size_t jj=0; jj<it->size(); ++jj)
//                 {
//                     if ((*it)[jj]<0)
//                         break;
//                     std::cout << (*it)[jj] << ",";
//                 }
//                 std::cout << ": " << curr_inters << ", " << area1 << ", " << area2 << std::endl;

                found = true;
                break;
            }
//             std::cout << " " << curr_inters << ", " << area1 << ", " << area2 << std::endl;
        }
        
        if (!found)
        {
            std::vector<int> to_put;
            for (std::size_t kk=0; kk<pairs[ii].size(); ++kk)
                to_put.push_back(pairs[ii][kk]);

            all_pairs.push_back(to_put);
        }
        
        
//         std::vector<unsigned int> to_put = leaves[pairs[ii][0]];
//         for (std::size_t jj=1; jj<n_neighs; ++jj)
//         {
//             if (pairs[ii][jj]<0)
//                 break;
//             to_put.insert(to_put.begin(),leaves[pairs[ii][jj]].begin(),leaves[pairs[ii][jj]].end());
//         }
//         std::sort(to_put.begin(), to_put.end());
//         all_pairs.insert(to_put);
        
//                 for (std::size_t jj=0; jj<to_put.size(); ++jj)
//         {
//             std::cout << to_put[jj]+1 << ", ";
//         }
//         std::cout << std::endl;
    }
    
    /* Assign output parameters */
    std::size_t n_pairs_2 = all_pairs.size();
    plhs[0] = mxCreateDoubleMatrix(n_pairs_2,n_neighs,mxREAL);
    MatlabMultiArray<double> pairs2(plhs[0]);
    
    // Fill pairs2
    std::size_t curr_pair = 0;
    std::list<std::vector<int> >::iterator list_it = all_pairs.begin();
    for (; list_it!=all_pairs.end(); ++list_it)
    {
        std::vector<int>::const_iterator vec_it2 = list_it->begin();
        std::size_t curr_reg = 0;
        for (; vec_it2!=list_it->end(); ++vec_it2)
        {
            pairs2[curr_pair][curr_reg] = *vec_it2+1;
            curr_reg++;
        }
        curr_pair++;
    }
    
}
