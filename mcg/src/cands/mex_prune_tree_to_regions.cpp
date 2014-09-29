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

#include <cstring>
#include <string>
#include <iostream>
#include <set>
#include <map>
#include <list>
#include <queue>
#include "matlab_multiarray.hpp"

/* Input Arguments */
#define	PART      prhs[0]
#define	MER_SEQ   prhs[1]
#define	REGS      prhs[2]

/* Output Arguments */
#define	NEW_PART    plhs[0]
#define	NEW_MER_SEQ plhs[1]
#define	LUT         plhs[2]

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
    ConstMatlabMultiArray<double> part(PART);
    ConstMatlabMultiArray<double> mer_seq(MER_SEQ);
    ConstMatlabMultiArray<double> regions(REGS);
        
    /* Input sizes and checks */
    std::size_t sm = part.shape()[0];
    std::size_t sn = part.shape()[1];
    std::size_t n_sons_max= mer_seq.shape()[1]-1;
    std::size_t n_merges = mer_seq.shape()[0];
    std::size_t n_leaves  = mer_seq[0][n_sons_max]; // n_merges+1;  --> Not valid for non-binary trees
    std::size_t n_regs    = n_leaves+n_merges;

    /* Allocate out partition */
    NEW_PART  =  mxCreateDoubleMatrix(sm, sn, mxREAL);
    MatlabMultiArray<double> new_part(NEW_PART);
        
    // Start descendants
    std::vector<std::vector<double> > descendants(n_regs);
    for (std::size_t ii=0; ii<n_regs; ++ii)
        descendants[ii].push_back(ii);
    
    /* Needed regions to set */
    std::set<double> needed_regions;
    if (regions.shape()[0]==1)
    {
        for(std::size_t ii=0; ii<regions.shape()[1]; ++ii)
        {
            needed_regions.insert(regions[0][ii]);
        }
    }
    else if (regions.shape()[1]==1)
    {
        for(std::size_t ii=0; ii<regions.shape()[0]; ++ii)
        {
            needed_regions.insert(regions[ii][0]);
        }
    }
    else
        mexErrMsgTxt("ERROR");
        
    /* Perform the maximum number of mergings that keep regions */
    std::list<std::size_t> mergs_to_keep;
    std::set<std::size_t> leaves_to_keep;
    for(std::size_t ii=0; ii<n_merges; ++ii)
    {
        double parent = mer_seq[ii][n_sons_max];
        bool merging_can_be_done = true;
        
        // Evolve merging sequence
        for (std::size_t jj=0; jj<n_sons_max; ++jj)
        {
            double curr_son = mer_seq[ii][jj];
            if (curr_son<0)
                break;

            // Gather descendants
            descendants[parent].insert(descendants[parent].begin(),descendants[curr_son].begin(),descendants[curr_son].end());
        
            // Check if current sons are needed regions
            if(needed_regions.find(curr_son)!=needed_regions.end())
            {
                merging_can_be_done = false;
                needed_regions.insert(parent);
                mergs_to_keep.push_back(ii);
                break;
            }
        }
        
        if (merging_can_be_done)
        {
            leaves_to_keep.insert(parent);
            for (std::size_t jj=0; jj<n_sons_max; ++jj)
            {
                double curr_son = mer_seq[ii][jj];
                if (curr_son<0)
                    break;
                leaves_to_keep.erase(curr_son);
            }
        }
    }
    
    /* Create pruned merging sequence */
    boost::array<std::size_t,2> dims; dims[0] = mergs_to_keep.size(); dims[1] = n_sons_max+1;
    boost::multi_array<double,2> curr_ms(dims);
    std::list<std::size_t>::iterator it=mergs_to_keep.begin();
    std::size_t curr_n_merges = 0;
    for(; it!=mergs_to_keep.end(); ++it)
    {
        for (std::size_t jj=0; jj<n_sons_max+1; ++jj)
        {
            curr_ms[curr_n_merges][jj] = mer_seq[*it][jj];
        }
        curr_n_merges++;
    }
      
            
    /* Create LUT*/
    typedef std::map<std::size_t,std::size_t> map_type;
    map_type lut;
    std::set<std::size_t>::iterator it2 = leaves_to_keep.begin();
    for( ; it2!=leaves_to_keep.end(); ++it2)
    {
        std::vector<double>::const_iterator it3 = descendants[*it2].begin();
        for(; it3!=descendants[*it2].end(); ++it3)
        {
            lut.insert(std::pair<std::size_t,std::size_t>(*it3,*it2));
        }
    }
    
    /* Perform unneeded merging on partition */
    dims[0] = sm; dims[1] = sn;
    boost::multi_array<double,2> curr_part(dims);
    for(std::size_t ii=0; ii<sm; ++ii)
    {
        for(std::size_t jj=0; jj<sn; ++jj)
        {
            map_type::const_iterator it4 = lut.find(part[ii][jj]);
            if(it4==lut.end())
            {
                curr_part[ii][jj] = part[ii][jj];
            }
            else
            {
                curr_part[ii][jj] = it4->second;
            }
        }
    }

    /* Relabel part and create LUT */
    typedef std::map<std::size_t,std::size_t> map_type;
    map_type curr_lut;
    std::size_t curr_reg = 0;
    for(std::size_t ii=0; ii<sm; ++ii)
    {
        for(std::size_t jj=0; jj<sn; ++jj)
        {
            map_type::iterator map_it = curr_lut.find(curr_part[ii][jj]);
            if (map_it==curr_lut.end())
            {
                curr_lut.insert(std::pair<std::size_t,std::size_t>(curr_part[ii][jj], curr_reg));
                new_part[ii][jj] = curr_reg + 1;
                curr_reg++;
            }
            else
            {
                new_part[ii][jj] = map_it->second + 1;
            }
        }
    }
    
    /* Relabel merging sequence */
    NEW_MER_SEQ =  mxCreateDoubleMatrix(curr_ms.shape()[0], curr_ms.shape()[1], mxREAL);
    MatlabMultiArray<double> new_ms(NEW_MER_SEQ);
    for(std::size_t ii=0; ii<new_ms.shape()[0]; ++ii)
    {
        double parent = curr_ms[ii][n_sons_max];
        curr_lut.insert(std::pair<std::size_t,std::size_t>(parent, curr_reg));
        new_ms[ii][n_sons_max] = curr_reg + 1;
        curr_reg++;
        
        // Evolve merging sequence
        for (std::size_t jj=0; jj<n_sons_max; ++jj)
        {
            double curr_son = curr_ms[ii][jj];
            if (curr_son<0)
                break;
            map_type::iterator map_it = curr_lut.find(curr_son);
            if (map_it==curr_lut.end())
            {
                mexErrMsgTxt("ERROR 2");
            }
            else
            {
                new_ms[ii][jj] = map_it->second+1;
            }
        }
    }
    
    /* Inverse LUT and sort it */
    map_type sort_lut;
    map_type::iterator it5 = curr_lut.begin();
    for(; it5!=curr_lut.end(); ++it5)
        sort_lut.insert(std::pair<std::size_t,std::size_t>(it5->second, it5->first));

    /* Store it */
    LUT =  mxCreateDoubleMatrix(sort_lut.size(), 2, mxREAL);
    MatlabMultiArray<double> out_lut(LUT);
    it5 = sort_lut.begin();
    for(std::size_t ii=0; it5!=sort_lut.end(); ++it5,++ii)
    {
        out_lut[ii][0] = it5->first + 1;
        out_lut[ii][1] = it5->second + 1;
    }
}

