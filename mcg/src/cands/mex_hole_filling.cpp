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

void vector_difference(std::vector<unsigned int>& vec1, std::vector<unsigned int>& vec2, std::vector<unsigned int>& vec_diff)
{
    vec_diff.resize(vec1.size());
    std::vector<unsigned int>::iterator it = std::set_difference(
            vec1.begin(), vec1.end(),
            vec2.begin(), vec2.end(), vec_diff.begin());
    vec_diff.resize(it-vec_diff.begin());
}

void vector_union(std::vector<unsigned int>& vec1, std::vector<unsigned int>& vec2, std::vector<unsigned int>& vec_union)
{
    vec_union.resize(vec1.size()+vec2.size());
    std::vector<unsigned int>::iterator it = std::set_union(
            vec1.begin(), vec1.end(),
            vec2.begin(), vec2.end(), vec_union.begin());
    vec_union.resize(it-vec_union.begin());
}

using namespace std;

void mexFunction( int nlhs, mxArray *plhs[], 
        		  int nrhs, const mxArray*prhs[] )
{    
    /* Check for proper number of arguments */
    if (nrhs != 5) { 
    	mexErrMsgTxt("Five input arguments required."); 
    } else if (nlhs > 5) {
        mexErrMsgTxt("Too many output arguments."); 
    } 
    
    /* Input as a Multiarray */
    ConstMatlabMultiArray<double> lp(prhs[0]); /* Leaves partition */
    ConstMatlabMultiArray<double> ms(prhs[1]); /* Merging sequence */
    ConstMatlabMultiArray<double> cands(prhs[2]); /* Candidates */
    ConstMatlabMultiArray<double> neigh_pairs_min(prhs[3]); 
    ConstMatlabMultiArray<double> neigh_pairs_max(prhs[4]);
    
    /* Input sizes and checks */
    size_t sm = lp.shape()[0];
    size_t sn = lp.shape()[1];
    size_t n_sons_max= ms.shape()[1]-1;
    size_t n_merges  = ms.shape()[0];
    size_t n_leaves  = ms[0][n_sons_max]-1; // n_merges+1;  --> Not valid for non-binary trees
    size_t n_regs    = n_leaves+n_merges;
    size_t n_cands   = cands.shape()[0];
    size_t n_regs_max= cands.shape()[1];
    size_t n_pairs   = neigh_pairs_min.shape()[0];

    
    /* Create LUT between region indices and merging  *
     * sequence positions where they disappear        */
    map<size_t,size_t> lut;
    for(size_t ii=0; ii<n_merges; ++ii)
    {
        for(size_t jj=0; jj<n_sons_max; ++jj)
        {
            if (ms[ii][jj]==0)
                break;
            lut[ms[ii][jj]] = ii;
        }
    }
    

    // *********************************************************************************
    //   Get all neighbors from each region at any levels of the hierarchy
    // *********************************************************************************
    // Prepare cells to store results
    vector<vector<unsigned int> >   neighbors(n_regs+1); // +1 to include border
    vector<vector<unsigned int> > descendants(n_regs+1);
    
    // Start leaves
    for (size_t ii=1; ii<=n_leaves; ++ii)
        descendants[ii].push_back(ii);

    // Add initial pairs
    for (size_t ii=0; ii<n_pairs; ++ii)
    {
        unsigned int min_id = neigh_pairs_min[ii][0];
        unsigned int max_id = neigh_pairs_max[ii][0];        
        neighbors[min_id].push_back(max_id);
        neighbors[max_id].push_back(min_id);
    }
    
    // Sort the neighbors (for efficiency later)
    for (size_t ii=0; ii<=n_regs; ++ii)
        sort(neighbors[ii].begin(),neighbors[ii].end());
    
    // Evolve through merging sequence
    for (size_t ii=0; ii<n_merges; ++ii)
    {
        unsigned int parent = ms[ii][n_sons_max];
        
        vector<unsigned int> all_neighs;
        vector<unsigned int> tmp_neighs;
        for (size_t jj=0; jj<n_sons_max; ++jj)
        {
            if (ms[ii][jj]==0)
                break;
            vector_union(all_neighs,neighbors[ms[ii][jj]],tmp_neighs);
            all_neighs = tmp_neighs;
        }
        
        // Descendants (including itself)
        descendants[parent].push_back(parent);
        for (size_t jj=0; jj<n_sons_max; ++jj)
        {
            if (ms[ii][jj]==0)
                break;
            descendants[parent].insert(descendants[parent].begin(),
                                       descendants[ms[ii][jj]].begin(),
                                       descendants[ms[ii][jj]].end());
        }
        sort(descendants[parent].begin(),descendants[parent].end());       
        vector_difference(all_neighs,descendants[parent],neighbors[parent]);
                
        // Update all neighbors
        for (size_t jj=0; jj<neighbors[parent].size(); ++jj)
        {
            unsigned int curr_neigh = neighbors[parent][jj];
            neighbors[curr_neigh].push_back(parent);
        }
    }
    // *********************************************************************

    
    
    // *********************************************************************
    //  Sweep and process all candidates - This could be parallelized
    // *********************************************************************
    size_t n_max_cands_hf_out   = n_regs_max;
    size_t n_max_cands_comp_out = 1;
    vector<list<double> > all_cands_hf(n_cands);
    vector<list<double> > all_cands_comp(n_cands);
    for(size_t ii=0; ii<n_cands; ++ii)
    {
        /* Is it the whole image? */
        if (cands[ii][0]==n_regs)
        {
            all_cands_hf[ii].push_back(n_regs);
            all_cands_comp[ii].push_back(n_regs); 
            // Complement of all image, it should be none, but this would give
            // us problems, we set all image, which will be removed in deduplication
        }
        else
        {
            /* -------------------------------------------------------*
             * Get the minimum set of regions needed for the regions  *
             *   in the candidates to coexist in the same partition,  *
             *   that is, 'coex' is the smallest set of labels so     *
             *   that the partition formed by these labels contains   *
             *   all regions in the current candidate                 *
             *--------------------------------------------------------*/
            set<double> cand_regs;// Regions in the current candidate
            set<double> coex;     // Coexistent regions
            set<double> coex_anc; // Ancestors of the coexistent regions

            // Get the minimum merging sequence index to start sweeping
            //  and fill the sets with the initial values
            size_t min_ms_index = n_merges+1;
            for(size_t jj=0; jj<n_regs_max; ++jj)
            {
                if (cands[ii][jj]==0)
                    break;
                min_ms_index = min(min_ms_index,lut[cands[ii][jj]]);

                cand_regs.insert(cands[ii][jj]);
                coex.insert(cands[ii][jj]);
                coex_anc.insert(cands[ii][jj]);
            }

            // Sweep
            for(size_t jj=min_ms_index; jj<n_merges; ++jj)
            {
                // Is the current merging needed?
                bool needed_merg = false;
                for(size_t kk=0; kk<n_sons_max; ++kk)
                {
                    if (ms[jj][kk]==0)
                        break;
                    if (coex_anc.count(ms[jj][kk]))
                        needed_merg = true;
                }

                // If needed add each region to the corresponding set
                if (needed_merg)
                {
                    // Parent to ancestors
                    coex_anc.insert(ms[jj][n_sons_max]);

                    // Siblings to coexistent
                    for(size_t kk=0; kk<n_sons_max; ++kk)
                    {
                        if (ms[jj][kk]==0)
                            break;
                        if (coex_anc.count(ms[jj][kk])==0)
                            coex.insert(ms[jj][kk]);
                    }
                }
            }
            /* -------------------------------------------------------*/

            
            
            // Remove the labels in the candidate from coexistent
            //  and add the border (0)
            coex.insert(0);
            for(size_t jj=0; jj<n_regs_max; ++jj)
            {
                if (cands[ii][jj]==0)
                    break;
                coex.erase(cands[ii][jj]);
            }

            
            /* -------------------------------------------------------*
             * Now form the sets of connected regions: those without  *
             *  connection to the border will be holes                *
             *--------------------------------------------------------*/
            vector<set<double> > curr_groups;
            set<double>::iterator it = coex.begin();
            for (; it!=coex.end(); ++it)
            {
                int found_group = -1;
                vector<size_t> to_erase;
                for(size_t jj=0; jj<curr_groups.size(); ++jj)
                {
                    // See if we can connect some neighbor of the current
                    //  coexistent region to any formed set
                    for(size_t kk=0; kk<neighbors[*it].size(); ++kk)
                    {
                        if (curr_groups[jj].count(neighbors[*it][kk]))
                        {
                            if (found_group<0)
                            {
                                // Add current coexistent to group
                                curr_groups[jj].insert(*it);
                            }
                            else
                            {
                                // Merge two groups linked by current coexistent
                                curr_groups[jj].insert(curr_groups[found_group].begin(),curr_groups[found_group].end());
                                to_erase.push_back(found_group);
                            }
                            found_group = jj;
                            break;
                        }
                    }
                }
                
                if (found_group<0)
                {
                    // Create new group
                    set<double> to_add;
                    to_add.insert(*it);
                    curr_groups.push_back(to_add);
                }
                else
                {
                    // Erase unneeded groups (reverse order to keep indices correct)
                    for(size_t jj=to_erase.size(); jj>0; --jj)
                    {
                        curr_groups.erase(curr_groups.begin()+to_erase[jj-1]);
                    }
                }
            }
             
            // Add current labels to new hole filling candidate
            size_t curr_n_regs_hf   = 0;
            for(size_t jj=0; jj<n_regs_max; ++jj)
            {
                if (cands[ii][jj]==0)
                    break;
                all_cands_hf[ii].push_back(cands[ii][jj]);
                curr_n_regs_hf++;
            }
            
            // Groups without 'touching' the border --> Holes
            // Groups 'touching' the border --> Complement
            size_t curr_n_regs_comp = 0;
            for(size_t jj=0; jj<curr_groups.size(); ++jj)
            {
                if (curr_groups[jj].count(0)) // Border, so add to complementary
                {
                    set<double>::iterator it = curr_groups[jj].begin();
                    for (; it!=curr_groups[jj].end(); ++it)
                    {
                        all_cands_comp[ii].push_back(*it);
                        curr_n_regs_comp++;
                    }                    
                }
                else // No border, so add to holes
                {
                    set<double>::iterator it = curr_groups[jj].begin();
                    for (; it!=curr_groups[jj].end(); ++it)
                    {
                        all_cands_hf[ii].push_back(*it);
                        curr_n_regs_hf++;
                    }
                }
            }
            n_max_cands_hf_out   = max(n_max_cands_hf_out,   curr_n_regs_hf);
            n_max_cands_comp_out = max(n_max_cands_comp_out, curr_n_regs_comp);
            /* -------------------------------------------------------*/
        }
    }
    
    /* Allocate out cands */
    plhs[0] = mxCreateDoubleMatrix(n_cands, n_max_cands_hf_out, mxREAL);
    MatlabMultiArray<double> cands_hf(plhs[0]);
    for(size_t ii=0; ii<n_cands; ++ii)
    {
        list<double>::iterator it = all_cands_hf[ii].begin();
        for (size_t jj=0; it!=all_cands_hf[ii].end(); ++it, ++jj)
        {
            cands_hf[ii][jj] = *it;
        }
    }
    
    if (n_max_cands_comp_out>1)
        n_max_cands_comp_out--; // -1 to ignore 0
                
    plhs[1] = mxCreateDoubleMatrix(n_cands, n_max_cands_comp_out, mxREAL); 
    MatlabMultiArray<double> cands_comp(plhs[1]);
    for(size_t ii=0; ii<n_cands; ++ii)
    {
        list<double>::iterator it = all_cands_comp[ii].begin();
        ++it; // Ignore 0
        if(it==all_cands_comp[ii].end())
        {
            cands_comp[ii][0] = n_regs; // If the complementary is null,
                                        // we set all the image (will be ignored when deduplicating)
        }
        else
        {
            for (size_t jj=0; it!=all_cands_comp[ii].end(); ++it, ++jj)
            {
                cands_comp[ii][jj] = *it;
            }
        }
    }
}

