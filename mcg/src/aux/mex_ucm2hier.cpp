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
#include <map>
#include <algorithm>

using namespace std;


struct cont_elem
{
    cont_elem(size_t new_x, size_t new_y, size_t new_n1_x, size_t new_n1_y, size_t new_n2_x, size_t new_n2_y)
                 : x(new_x), y(new_y), n1_x(new_n1_x), n1_y(new_n1_y), n2_x(new_n2_x), n2_y(new_n2_y)
    {
    }

    size_t x;
    size_t y;
    size_t n1_x;
    size_t n1_y;
    size_t n2_x;
    size_t n2_y;
};

void insert_labels_to_merge(double lab1,double lab2,list<set<double> >& to_merge)
{
    list<set<double> >::iterator l1 = to_merge.begin();
    for( ; l1!=to_merge.end(); ++l1)
        if (l1->find(lab1) != l1->end())
            break;
    list<set<double> >::iterator l2 = to_merge.begin();
    for( ; l2!=to_merge.end(); ++l2)
        if (l2->find(lab2) != l2->end())
            break;

    if (l1==to_merge.end() && l2==to_merge.end())  // Neither label found --> New set
    {
        set<double> to_put;
         to_put.insert(lab1);
         to_put.insert(lab2);
        to_merge.push_back(to_put);
    }
    else if (l1==to_merge.end() && l2!=to_merge.end())  // Found one of the two in one set --> Add the other
        l2->insert(lab1);
    else if (l1!=to_merge.end() && l2==to_merge.end())  // Found one of the two in one set --> Add the other
        l1->insert(lab2);
    else if (l1!=l2)  // Both found in different sets --> Merge the two sets
    {
        l1->insert(l2->begin(), l2->end());
        to_merge.erase(l2);
    }
}


void mexFunction( int nlhs, mxArray *plhs[], 
        		  int nrhs, const mxArray*prhs[] )
{
    if(nrhs!=2)
        mexErrMsgTxt("There should be 2 input parameters");
    
    /* Input parameters */
    ConstMatlabMultiArray<double> lp(prhs[0]);
    ConstMatlabMultiArray<double> ucm2(prhs[1]);

    std::size_t sx   = lp.shape()[0];
    std::size_t sy   = lp.shape()[1];

    
    /*-------------------------------------------------------------*/
    /* Create LUT of contour positions and sets of labels to merge */
    /*-------------------------------------------------------------*/
    typedef map<double,list<set<double> > > map_type;
    map_type ucm_th_leaves_pairs;
    
    /* Scan horizontal contours*/
    for (std::size_t xx=2; xx<2*sx; xx+=2)
    {
        for (std::size_t yy=1; yy<2*sy; yy+=2)
        {
            if (ucm2[xx][yy]>0)
            {
                cont_elem cont(xx,yy,(xx/2)-1,(yy-1)/2,(xx/2),(yy-1)/2);
                double lab1 = min(lp[cont.n1_x][cont.n1_y],lp[cont.n2_x][cont.n2_y]);
                double lab2 = max(lp[cont.n1_x][cont.n1_y],lp[cont.n2_x][cont.n2_y]);

                /* Get the sets of regions to be merged                       */
                /* In one threshold there can be more than one region forming */
                map_type::iterator map_it = ucm_th_leaves_pairs.find(ucm2[xx][yy]);
                if (map_it==ucm_th_leaves_pairs.end()) // New ucm threshold
                {
                    set<double> set_to_put;
                     set_to_put.insert(lab1);
                     set_to_put.insert(lab2);
                    list<set<double> > to_put;
                    to_put.push_back(set_to_put); // Add up and down neighbors
                    ucm_th_leaves_pairs.insert(map_type::value_type(ucm2[xx][yy], to_put));
                }
                else
                    insert_labels_to_merge(lab1,lab2,map_it->second);
            }
        }
    }

    /* Scan vertical contours*/
    for (std::size_t xx=1; xx<2*sx; xx+=2)
    {
        for (std::size_t yy=2; yy<2*sy; yy+=2)
        {
            if (ucm2[xx][yy]>0)
            {
                cont_elem cont(xx,yy,(xx-1)/2,(yy/2)-1,(xx-1)/2,(yy/2));

                double lab1 = min(lp[cont.n1_x][cont.n1_y],lp[cont.n2_x][cont.n2_y]);
                double lab2 = max(lp[cont.n1_x][cont.n1_y],lp[cont.n2_x][cont.n2_y]);

                /* Get the sets of regions to be merged                       */
                /* In one threshold there can be more than one region forming */
                map_type::iterator map_it = ucm_th_leaves_pairs.find(ucm2[xx][yy]);
                if (map_it==ucm_th_leaves_pairs.end()) // New ucm threshold
                {
                    set<double> set_to_put;
                     set_to_put.insert(lab1);
                     set_to_put.insert(lab2);
                    list<set<double> > to_put;
                    to_put.push_back(set_to_put); // Add left and right neighbors
                    ucm_th_leaves_pairs.insert(map_type::value_type(ucm2[xx][yy], to_put));
                }
                else
                    insert_labels_to_merge(lab1,lab2,map_it->second);
            }
        }
    }
    
    /* Number of thresholds is number of mergings */
    int n_merges = ucm_th_leaves_pairs.size();
    /*-------------------------------------------------------------*/
    
   
    /*-----------------------------------------------------------*/
    /*  Scan all thresholds of the UCM in increasing order and   */
    /*   do the mergings of the N regions that 'disappear'       */
    /*-----------------------------------------------------------*/

    /* Get number of leaves */
    double curr_max = 0;
    for (std::size_t xx=0; xx<sx; ++xx)
        for (std::size_t yy=0; yy<sy; ++yy)
            curr_max = max(curr_max, lp[xx][yy]);
    double n_leaves = curr_max;
    
    /* Start LUT with leaves */
    vector<double> lut(n_leaves);
    vector<set<double> > ilut(n_leaves);
    for (size_t ll=0; ll<n_leaves; ++ll)
    {
        lut[ll] = ll+1;
        ilut[ll].insert(ll+1);
    }
    
    /* Scan */
    size_t n_max_children = 1;
    vector<double> parent_labels;
    vector<list<double> > children_labels;
    vector<double> start_ths(n_leaves,0);
    
    map_type::iterator it = ucm_th_leaves_pairs.begin();
    for ( ; it!=ucm_th_leaves_pairs.end(); ++it)
    {        
        list<set<double> >& leaves_to_merge = it->second;
        list<set<double> >    regs_to_merge;
        
        list<set<double> >::iterator it3 = leaves_to_merge.begin();
        for ( ; it3!=leaves_to_merge.end(); ++it3)
        {
            /* Get the labels of the children regions from the leave labels */
            set<double>::iterator it4  = it3->begin();
            set<double>::iterator it4b = it3->begin();
            ++it4b;
            
            for( ; it4b!=it3->end(); ++it4b)
                insert_labels_to_merge(lut[*it4-1],lut[*it4b-1],regs_to_merge);
        }
        
        it3 = regs_to_merge.begin();
        for ( ; it3!=regs_to_merge.end(); ++it3)
        {
            curr_max += 1;
            ilut.push_back(set<double>());

            /* Update LUT */
            set<double> to_update(*it3);
            while(!to_update.empty())
            {
                /* Get the first to update */
                double curr = *(to_update.begin());             
                to_update.erase(to_update.begin());
                if (lut[curr-1]!=curr)
                    to_update.insert(lut[curr-1]);
                
                set<double>::iterator it5 = ilut[curr-1].begin();
                for( ; it5!=ilut[curr-1].end(); ++it5)
                {
                    lut[*it5-1] = curr_max;
                    ilut[curr_max-1].insert(*it5);
                }
                ilut[curr-1].clear();
            }

            /* Copy the labels to the merging sequence */
            children_labels.push_back(list<double>(it3->begin(),it3->end()));
            n_max_children = max(n_max_children,it3->size());
            
            /* Add parent to the LUT and to out*/
            lut.push_back(curr_max);
            parent_labels.push_back(curr_max);
            
            /* Store the start contour threshold */
            start_ths.push_back(it->first);
        }
    }
    double n_regs = curr_max;
    /*-----------------------------------------------------------*/
    
    if(n_regs!=start_ths.size())
        mexErrMsgTxt("Oh oh");
   
    /* Output allocation */
    plhs[0] = mxCreateDoubleMatrix(parent_labels.size(),n_max_children+1,mxREAL);
    MatlabMultiArray<double> ms_out(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(n_regs,1,mxREAL);
    MatlabMultiArray<double> start_ths_out(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(n_regs,1,mxREAL);
    MatlabMultiArray<double> end_ths_out(plhs[2]);
    
    /* Copy data to output */
    for (size_t ii=0; ii<n_regs; ++ii)
        start_ths_out[ii][0] = start_ths[ii];
        
    for (size_t ii=0; ii<parent_labels.size(); ++ii)
    {
        ms_out[ii][n_max_children] = parent_labels[ii];
        list<double>::iterator it4 = children_labels[ii].begin();
        for (size_t jj=0 ; it4!=children_labels[ii].end(); ++it4, ++jj)
        {
            ms_out[ii][jj] = *it4;
            end_ths_out[*it4-1][0] = start_ths_out[parent_labels[ii]-1][0];
        }
    }
    end_ths_out[n_regs-1][0] = 1;
}
    
