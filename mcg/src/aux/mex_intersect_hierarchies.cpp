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
#include <map>
#include <set>
#include <iostream>
#include "boost/multi_array.hpp"

/* Arguments */
#define	PARTS	 prhs[0]
#define	MERSEQS	 prhs[1]
#define	INT_PART plhs[0]
#define	LUT      plhs[1]

void mexFunction( int nlhs, mxArray *plhs[], 
        		  int nrhs, const mxArray*prhs[] )
{  
    int num_dims = mxGetNumberOfDimensions(PARTS);
    const int* in_dims = mxGetDimensions(PARTS);
    
    std::size_t num_parts;
    if (num_dims==3)
    	num_parts = in_dims[2];
    else if (num_dims==2)
        num_parts = 1;
    else
        mexErrMsgTxt("Input partitions must be 3- or 2-dimensional");
    
    std::size_t s_x = in_dims[0];
    std::size_t s_y = in_dims[1];

    if(mxGetClassID(PARTS)!=mxUINT32_CLASS)
        mexErrMsgTxt("Input must be UINT32");
    
    /* Input as a Multiarray */
    boost::multi_array_ref<int, 3> parts((int*)mxGetPr(PARTS), boost::extents[in_dims[0]][in_dims[1]][num_parts], boost::fortran_storage_order());
    
    /* Handy typedefs */
    typedef std::vector<int> vector_type;
    typedef std::set<int> set_type;
    typedef std::map<vector_type, int> map_type;
    typedef std::pair<vector_type, int> pair_type;
    typedef std::map<int,set_type> imap_type;
    typedef std::pair<int,set_type> ipair_type;
        
    /* Create output intersection partition */
    INT_PART = mxCreateNumericMatrix(s_x, s_y, mxUINT32_CLASS, mxREAL);
    boost::multi_array_ref<int, 2> int_part((int*)mxGetPr(INT_PART), boost::extents[s_x][s_y], boost::fortran_storage_order());
    
    /* Scan all partitions */
    int num_regions;
    map_type out_lut;
    imap_type label_lut;
    int next_label = 1;
    int curr_label = 0;
    for(std::size_t xx=0; xx<s_x; ++xx)
    {
        for(std::size_t yy=0; yy<s_y; ++yy)
        {
            /* Build current id vector */
            vector_type curr_ids(num_parts);
            for(std::size_t pp=0; pp<num_parts; ++pp)
            {
                curr_ids[pp] = parts[xx][yy][pp];
            }
            
            /* Have we already found it? */
            map_type::const_iterator it = out_lut.find(curr_ids);
            if (it==out_lut.end())
            {   
                curr_label = next_label;
                next_label++;
               
                out_lut.insert(pair_type(curr_ids, curr_label));       
            }
            else
            {
                curr_label = (*it).second;
            }
            int_part[xx][yy] = curr_label;
            
            /* Add to label lut */
            for(std::size_t pp=0; pp<num_parts; ++pp)
            {
                imap_type::iterator iit = label_lut.find(curr_ids[pp]);
                if (iit==label_lut.end())
                {
                    set_type to_insert;
                    to_insert.insert(curr_label);
                    label_lut.insert(ipair_type(curr_ids[pp],to_insert));
                }
                else
                {
                    (*iit).second.insert(curr_label);
                }
            }
        }
    }
    num_regions = out_lut.size();

    /* Invert label LUT */
    std::map<int, vector_type> ilabel_lut;
    for(map_type::const_iterator it=out_lut.begin(); it!=out_lut.end(); ++it)
    {
        ilabel_lut.insert(std::pair<int, vector_type>((*it).second, (*it).first));
    }
    
    /* Perform the mergings, if provided */
    if(nrhs==2)
    {
        if (!mxIsCell(MERSEQS))
            mexErrMsgTxt("Second argument mer_seqs must be a cell");

        int num_hierarchies = mxGetNumberOfElements(MERSEQS);
        
        for(std::size_t ii=0; ii<num_hierarchies; ++ii)
        {
            const mxArray* curr_mer_seq = mxGetCell(MERSEQS,ii);
            int num_mergings = mxGetNumberOfElements(curr_mer_seq);
            if (!mxIsStruct(curr_mer_seq))
                mexErrMsgTxt("Merging sequence bad formed: not a struct");
            for(std::size_t jj=0; jj<num_mergings; ++jj)
            {                  
                /* Get children */
                mxArray *children = mxGetField(curr_mer_seq, jj,"children");
                if (!mxIsDouble(children))
                    mexErrMsgTxt("Mergings must be doubles");               
                double* children_pr = mxGetPr(children);
                      
                /* Get parent */
                mxArray *parent = mxGetField(curr_mer_seq, jj,"parent");
                if (!mxIsDouble(parent))
                    mexErrMsgTxt("Parents must be doubles");               
                int parent_label = (int)*mxGetPr(parent);
  
                /* Create parent set for the lut */
                imap_type::iterator iit;
                set_type::iterator set_it;
                set_type parent_set;
                for(std::size_t kk=0; kk<mxGetN(children); ++kk)
                {
                    iit = label_lut.find((int)children_pr[kk]);
                    if (iit==label_lut.end())
                        mexErrMsgTxt("Bad formed tree: children not found at lut");
                    for(set_it = (*iit).second.begin(); set_it != (*iit).second.end(); ++set_it)
                        parent_set.insert(*set_it);
                }
                    
                /* Add parent entry at label_lut */
                iit = label_lut.find(parent_label);
                if (iit!=label_lut.end())
                    mexErrMsgTxt("Bad formed tree: parent already found at lut");
                label_lut.insert(ipair_type(parent_label,parent_set));
                
                /* Add parent label at each place marked by label_lut */
                for(set_it = parent_set.begin(); set_it != parent_set.end(); ++set_it)
                {
                    std::map<int, vector_type>::iterator it = ilabel_lut.find(*set_it);
                    if (it==ilabel_lut.end())
                        mexErrMsgTxt("The label should be at the lut");
                    (*it).second.push_back(parent_label);
                }
            }
        }
    }
    
    /* Store the LUT as a cell */
    int dims_out[1]; 
    dims_out[0] = num_regions;
    LUT = mxCreateCellArray(1, dims_out);
    
    for(std::map<int, vector_type>::const_iterator it=ilabel_lut.begin(); it!=ilabel_lut.end(); ++it)
    {
        mxArray* curr_entry = mxCreateNumericMatrix(1, (*it).second.size(), mxUINT32_CLASS, mxREAL);
        int* curr_pr = (int*)mxGetPr(curr_entry);
        std::size_t curr_jj = 0;
        for(vector_type::const_iterator vec_it=(*it).second.begin(); vec_it!=(*it).second.end(); ++vec_it)
        {
            curr_pr[curr_jj] = *vec_it;
            curr_jj++;
        }
        mxSetCell(LUT, (std::size_t)((*it).first-1), curr_entry);
    }
}
