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

void vector_difference(std::vector<unsigned int>& vec1, std::vector<unsigned int>& vec2, std::vector<unsigned int>& vec_diff)
{
    vec_diff.resize(vec1.size());
    std::vector<unsigned int>::iterator it = std::set_difference(
            vec1.begin(), vec1.end(),
            vec2.begin(), vec2.end(), vec_diff.begin());
    vec_diff.resize(it-vec_diff.begin());
}


void mexFunction( int nlhs, mxArray *plhs[], 
        		  int nrhs, const mxArray*prhs[] )
{   
    if(nrhs!=5)
        mexErrMsgTxt("There should be exactly 6 input parameters");
    
    /* Input parameters */
    ConstMatlabMultiArray<double> leaves_part(prhs[0]);
    ConstMatlabMultiArray<double> mer_seq    (prhs[1]);
    ConstMatlabMultiArray<double> neigh_pairs_min(prhs[2]);
    ConstMatlabMultiArray<double> neigh_pairs_max(prhs[3]);
    ConstMatlabMultiArray<double> num_cands_in(prhs[4]); // Number of pairs, triplets, etc.
    
    double n_reg_cand;
    std::vector<double> num_cands;
    if (num_cands_in.shape()[0]==1)
    {
        n_reg_cand = num_cands_in.shape()[1] + 1;
        num_cands.resize(n_reg_cand);
        for (std::size_t ii=1; ii<n_reg_cand; ++ii)
            num_cands[ii] = num_cands_in[0][ii-1];
    }
    else if (num_cands_in.shape()[1]==1)
    {
        n_reg_cand = num_cands_in.shape()[0] + 1;
        num_cands.resize(n_reg_cand);
        for (std::size_t ii=1; ii<n_reg_cand; ++ii)
            num_cands[ii] = num_cands_in[ii-1][0];
    }
    else
        mexErrMsgTxt("Number of candidates should be a vector");
     
    std::size_t sx         = leaves_part.shape()[0];
    std::size_t sy         = leaves_part.shape()[1];
    std::size_t n_merges   = mer_seq.shape()[0];
    std::size_t n_sons_max = mer_seq.shape()[1]-1;
    std::size_t n_leaves   = mer_seq[0][n_sons_max]; // n_merges+1;  --> Not valid for non-binary trees
    std::size_t n_regs     = n_leaves+n_merges;
    std::size_t n_pairs    = neigh_pairs_min.shape()[0];

    // *********************************************************************************
    //   Compute the neighbors, descendants, siblings, and ascendants of each region
    // *********************************************************************************
    // Prepare cells to store results
    std::vector<std::vector<unsigned int> > neighbors(n_regs);
    std::vector<std::vector<unsigned int> > descendants(n_regs);
    std::vector<std::vector<unsigned int> > ascendants(n_regs);
    std::vector<std::vector<unsigned int> > siblings(n_regs);
    
    // Start leaves
    for (std::size_t ii=0; ii<n_leaves; ++ii)
        descendants[ii].push_back(ii);

    // Add initial pairs
    for (std::size_t ii=0; ii<n_pairs; ++ii)
    {
        unsigned int min_id = neigh_pairs_min[ii][0];
        unsigned int max_id = neigh_pairs_max[ii][0];        
        neighbors[min_id].push_back(max_id);
        neighbors[max_id].push_back(min_id);
    }
    
    // Sort the neighbors (for efficiency later)
    for (std::size_t ii=0; ii<n_regs; ++ii)
        std::sort(neighbors[ii].begin(),neighbors[ii].end());
    
    // Evolve through merging sequence
    for (std::size_t ii=0; ii<n_merges; ++ii)
    {
        unsigned int parent = mer_seq[ii][n_sons_max];
        
        std::vector<unsigned int> all_neighs;
        std::vector<unsigned int> tmp_neighs;
        for (std::size_t jj=0; jj<n_sons_max; ++jj)
        {
            if (mer_seq[ii][jj]<0)
                break;
            vector_union(all_neighs,neighbors[mer_seq[ii][jj]],tmp_neighs);
            all_neighs = tmp_neighs;
        }
        
        // Descendants (including itself)
        descendants[parent].push_back(parent);
        for (std::size_t jj=0; jj<n_sons_max; ++jj)
        {
            if (mer_seq[ii][jj]<0)
                break;
            descendants[parent].insert(descendants[parent].begin(),
                                       descendants[mer_seq[ii][jj]].begin(),
                                       descendants[mer_seq[ii][jj]].end());
        }
        std::sort(descendants[parent].begin(),descendants[parent].end());       
        vector_difference(all_neighs,descendants[parent],neighbors[parent]);
                
        // Update all neighbors
        for (std::size_t jj=0; jj<neighbors[parent].size(); ++jj)
        {
            unsigned int curr_neigh = neighbors[parent][jj];
            neighbors[curr_neigh].push_back(parent);
        }
    }
        
    // Siblings
    for (std::size_t ii=0; ii<n_merges; ++ii)
    {
        // Get all regions that are merged at that step
        std::vector<unsigned int> all_siblings;
        for (std::size_t jj=0; jj<n_sons_max; ++jj)
        {
            if (mer_seq[ii][jj]<0)
                break;
            all_siblings.push_back(mer_seq[ii][jj]);
        }
        std::sort(all_siblings.begin(),all_siblings.end());       

        // For each region, add the rest of siblings
        for (std::size_t jj=0; jj<n_sons_max; ++jj)
        {
            if (mer_seq[ii][jj]<0)
                break;
            
            // siblings[mer_seq[ii][jj]] = all_siblings \ curr_reg 
            std::vector<unsigned int> curr_reg(1,mer_seq[ii][jj]);
            vector_difference(all_siblings,curr_reg,siblings[mer_seq[ii][jj]]);
        }
    }
    
    // Ascendants
    for (std::size_t ii=n_merges; ii>0; --ii)
    {
        unsigned int parent = mer_seq[ii-1][n_sons_max];
    
        for (std::size_t jj=0; jj<n_sons_max; ++jj)
        {
            if (mer_seq[ii-1][jj]<0)
                break;
            ascendants[mer_seq[ii-1][jj]].push_back(parent);
            
            std::vector<unsigned int> tmp;
            vector_union(ascendants[mer_seq[ii-1][jj]], ascendants[parent], tmp);
            ascendants[mer_seq[ii-1][jj]] = tmp;
        }
    }
    // *********************************************************************
    

    
    // *********************************************************************
    //           Top-down computation of pairs, triplets, etc.
    // *********************************************************************
    // All new singletons [0], pairs [1], triplets [2], etc.
    std::vector<std::list<std::vector<unsigned int> > > cands_list(n_reg_cand);
    std::vector<std:: set<std::vector<unsigned int> > > cands_set(n_reg_cand); // Set to remove duplicates
    std::vector<unsigned int> coexistent;
    unsigned int curr_n_reg_max = n_reg_cand;
    for (std::size_t ii=0; ii<n_merges; ++ii)
    {
        std::size_t curr_id = n_merges-ii-1;
                
        // Update coexistent
        //  1-Remove parent
        std::vector<unsigned int> tmp_coexistent;
        unsigned int parent = mer_seq[curr_id][n_sons_max];
        std::vector<unsigned int> tmp_parent(1,parent);
        vector_difference(coexistent,tmp_parent,tmp_coexistent);
        coexistent = tmp_coexistent;
        //  2-Add children
        for (std::size_t jj=0; jj<n_sons_max; ++jj)  
        {
            double child = mer_seq[curr_id][jj];
            if (child<0)
                break;
            coexistent.push_back(child);
        }
        std::sort(coexistent.begin(),coexistent.end());

        // All new singletons [0], pairs [1], and triplets [2]
        std::vector<std::list<std::vector<unsigned int> > > new_cands_list(n_reg_cand);
        std::vector<std:: set<std::vector<unsigned int> > > new_cands_set(n_reg_cand);

        // Add new singletons (children)
        for (std::size_t jj=0; jj<n_sons_max; ++jj)
        {
            double child = mer_seq[curr_id][jj];
            if (child<0)
                break;

            std::vector<unsigned int> to_put(1,child);
            new_cands_list[0].push_back(to_put);
            new_cands_set[0].insert(to_put);
            cands_list[0].push_back(to_put);
            cands_set[0].insert(to_put);
        }
        
        // Scan singletons to create pairs and pairs to create triplets (if needed), etc.
        for (std::size_t n_reg_id=1; n_reg_id<curr_n_reg_max; ++n_reg_id)
        {
            std::list<std::vector<unsigned int> >::iterator list_it = new_cands_list[n_reg_id-1].begin();
            for ( ; list_it!=new_cands_list[n_reg_id-1].end(); ++list_it)
            {
                // Regions forming the current candidate. We add regions to them to get from pairs to triplets or from singletons to pairs.
                std::vector<unsigned int> curr_regs_vec(*list_it);
                std::set   <unsigned int> curr_regs_set(curr_regs_vec.begin(),curr_regs_vec.end());
                
                /*******   Up_neighs:  All neighbors minus the descendants of all coexistent,       ********/
                /*******        to keep the order of the candidates in the hierarchy.               ********/
                /*******  We also remove the siblings, since they would create a repeated candidate ********/
                std::vector<unsigned int>  up_neighs;
                std::vector<unsigned int> tmp_neighs;

                // Add all neighbors from all regions (and then remove own and descendants)
                up_neighs = neighbors[curr_regs_vec[0]];
                for (std::size_t kk=1; kk<curr_regs_vec.size(); ++kk)
                {
                    vector_union(up_neighs,neighbors[curr_regs_vec[kk]],tmp_neighs);
                    up_neighs = tmp_neighs;
                }
                std::sort(up_neighs.begin(),up_neighs.end());
                // Remove own, descendants, and ascendants
                for (std::size_t kk=0; kk<curr_regs_vec.size(); ++kk)
                {
                    vector_difference( up_neighs,descendants[curr_regs_vec[kk]],tmp_neighs);
                    vector_difference(tmp_neighs, ascendants[curr_regs_vec[kk]], up_neighs);
                }

                // Scan all coexistent
                for (std::size_t kk=0; kk<coexistent.size(); ++kk)
                {
                    double coex = coexistent[kk];
                    if (curr_regs_set.find(coex)==curr_regs_set.end())
                    {
                        // Get descendants without own coexistent (to keep it as a neighbor)
                        // curr_descendants = descendants[coex] \ curr_coex
                        std::vector<unsigned int> curr_descendants;
                        std::vector<unsigned int> curr_coex(1,coex);
                        vector_difference(descendants[coex],curr_coex,curr_descendants);

                        // Remove descendants from current coexistent
                        // up_neighs = up_neighs \ curr_descendants
                        std::vector<unsigned int> tmp_neighs;
                        vector_difference(up_neighs,curr_descendants,tmp_neighs);
                        up_neighs = tmp_neighs;
                    }
                }
                
                // Remove siblings (just if they are pairs, otherwise we would be missing some of them)
                std::vector<unsigned int> curr_siblings;
                for (std::size_t kk=0; kk<curr_regs_vec.size(); ++kk)
                    if (siblings[curr_regs_vec[kk]].size()==1)
                        curr_siblings.push_back(siblings[curr_regs_vec[kk]][0]);
                std::sort(curr_siblings.begin(),curr_siblings.end());
                
                std::vector<unsigned int> tmp_neighs2;
                vector_difference(up_neighs,curr_siblings,tmp_neighs2);
                up_neighs = tmp_neighs2;
                
                /******************** Up_neighs updated ********************/
                
                
                // Store all up_neighs U current regions
                for (std::size_t kk=0; kk<up_neighs.size(); ++kk)
                {
                    // to_put = curr_regs_vec U up_neighs[up_neighs.size()-kk-1]
                    std::vector<unsigned int> to_put(curr_regs_vec);
                    to_put.push_back(up_neighs[up_neighs.size()-kk-1]);
                    
                    std::sort(to_put.begin(),to_put.end());

                    if (new_cands_set[n_reg_id].find(to_put)==new_cands_set[n_reg_id].end())
                    {
                        new_cands_list[n_reg_id].push_back(to_put);
                        new_cands_set [n_reg_id].insert(to_put);
                    }
                    if (cands_set[n_reg_id].find(to_put)==cands_set[n_reg_id].end())
                    {
                        cands_list[n_reg_id].push_back(to_put);
                        cands_set [n_reg_id].insert(to_put);
                    }
                }
            }
        }
        
        // Update curr_n_reg_max
        bool done = true;
        for (std::size_t n_reg_id=n_reg_cand-1; n_reg_id>0; --n_reg_id)
        {
            if (cands_set[n_reg_id].size()<num_cands[n_reg_id])
                done = false;
            else if (done)
                curr_n_reg_max = n_reg_id;
        }
        
        // Are we done?
        if (done)
            break;
    }
    // *********************************************************************

    // Store at output variable cell
    plhs[0]=mxCreateCellMatrix(n_reg_cand-1, 1);
    for (std::size_t kk=1; kk<n_reg_cand; ++kk)
    {       
        // Allocate the space at each slot of the cell
        std::size_t max_num_cands = cands_set[kk].size();
        double curr_num_cands = std::min((double)max_num_cands,(double)num_cands[kk]);
        mxArray *curr_cands = mxCreateDoubleMatrix(curr_num_cands,kk+1,mxREAL);
        MatlabMultiArray<double> cands_out(curr_cands);
        
        // Copy result to output
        std::list<std::vector<unsigned int> >::const_iterator list_it = cands_list[kk].begin();
        for (std::size_t ii=0; ii<curr_num_cands; ++ii)
        {
            for (std::size_t jj=0; jj<kk+1; ++jj)
            {
                cands_out[ii][jj] = (*list_it)[jj] + 1;
            }
            ++list_it;
        }
        
        // Set each slot of the cell
        mxSetCell(plhs[0],kk-1,curr_cands);
    }
}

