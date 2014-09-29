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


void mexFunction( int nlhs, mxArray *plhs[], 
        		  int nrhs, const mxArray*prhs[] )
{
    if(nrhs!=3)
        mexErrMsgTxt("There should be 3 input parameters");
    
    /* Input parameters */
    ConstMatlabMultiArray<double> lp(prhs[0]);
    ConstMatlabMultiArray<double> ucm2(prhs[1]);
    double nreg_max = mxGetScalar(prhs[2]);

    std::size_t sx   = lp.shape()[0];
    std::size_t sy   = lp.shape()[1];

    /* Output allocation */
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    MatlabMultiArray<double> ucm_th(plhs[0]);
    
    /*-----------------------------------------------------------*/
    /* Create LUT of contour positions where each ucm_th appears */
    /*-----------------------------------------------------------*/
    typedef map<double,list<cont_elem> > map_type;
    map_type cont_pos;
    
    /* Scan horizontal contours*/
    for (std::size_t xx=2; xx<2*sx; xx+=2)
    {
        for (std::size_t yy=1; yy<2*sy; yy+=2)
        {
            if (ucm2[xx][yy]>0)
            {
                map_type::iterator map_it = cont_pos.find(ucm2[xx][yy]);
                if (map_it==cont_pos.end()) // New label
                {
                    list<cont_elem> to_put;
                    to_put.push_back(cont_elem(xx,yy,(xx/2)-1,(yy-1)/2,(xx/2),(yy-1)/2)); // Add up and down neighbors
                    cont_pos.insert(map_type::value_type(ucm2[xx][yy], to_put));
                }
                else
                    map_it->second.push_back(cont_elem(xx,yy,(xx/2)-1,(yy-1)/2,(xx/2),(yy-1)/2));
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
                map_type::iterator map_it = cont_pos.find(ucm2[xx][yy]);
                if (map_it==cont_pos.end()) // New label
                {
                    list<cont_elem> to_put;
                    to_put.push_back(cont_elem(xx,yy,(xx-1)/2,(yy/2)-1,(xx-1)/2,(yy/2))); // Add up and down neighbors
                    cont_pos.insert(map_type::value_type(ucm2[xx][yy], to_put));
                }
                else
                    map_it->second.push_back(cont_elem(xx,yy,(xx-1)/2,(yy/2)-1,(xx-1)/2,(yy/2)));
            }
        }
    }
    
    /* Number of thresholds is number of mergings */
    int n_merges = cont_pos.size();
    /*-----------------------------------------------------------*/
    
    
   
    /*-----------------------------------------------------------*/
    /*  Scan all thresholds of the UCM in increasing order and   */
    /*   do the mergings of the N regions that 'disappear'       */
    /*-----------------------------------------------------------*/

    /* Copy lp and get number of regions */
    double curr_max = 0;
    boost::array<int,2> dims;  dims[0] = sx;  dims[1] = sy;
    boost::multi_array<double,2> curr_lp(dims);
    for (std::size_t xx=0; xx<sx; ++xx)
        for (std::size_t yy=0; yy<sy; ++yy)
        {
            curr_lp[xx][yy] = lp[xx][yy];
            curr_max = max(curr_max, lp[xx][yy]);
        }
    double curr_nl = curr_max;
    
    /* Scan */
    list<pair<double,double> > num_regs_at_th;
    num_regs_at_th.push_back(pair<double,double>(0,curr_nl));
    map_type::iterator it = cont_pos.begin();
    for ( ; it!=cont_pos.end(); ++it)
    {        
        /* Get the sets of regions to be merged                       */
        /* In one threshold there can be more than one region forming */
        list<set<double> > to_merge;
        list<cont_elem>::iterator it2=it->second.begin();
        for ( ; it2!=it->second.end(); ++it2)
        {   
            list<set<double> >::iterator l1 = to_merge.begin();
            for( ; l1!=to_merge.end(); ++l1)
                if (l1->find(curr_lp[it2->n1_x][it2->n1_y]) != l1->end())
                    break;
            list<set<double> >::iterator l2 = to_merge.begin();
            for( ; l2!=to_merge.end(); ++l2)
                if (l2->find(curr_lp[it2->n2_x][it2->n2_y]) != l2->end())
                    break;
            
            if (l1==to_merge.end() && l2==to_merge.end())  // Neither label found --> New set
            {
                set<double> to_put;
                to_put.insert(curr_lp[it2->n1_x][it2->n1_y]);
                to_put.insert(curr_lp[it2->n2_x][it2->n2_y]);
                to_merge.push_back(to_put);
            }
            else if (l1==to_merge.end() && l2!=to_merge.end())  // Found one of the two in one set --> Add the other
                l2->insert(curr_lp[it2->n1_x][it2->n1_y]);
            else if (l1!=to_merge.end() && l2==to_merge.end())  // Found one of the two in one set --> Add the other
                l1->insert(curr_lp[it2->n2_x][it2->n2_y]);
            else if (l1!=l2)  // Both found in different sets --> Merge the two sets
            {
                l1->insert(l2->begin(), l2->end());
                to_merge.erase(l2);
            }
        }
        
        list<set<double> >::iterator it3 = to_merge.begin();
        for ( ; it3!=to_merge.end(); ++it3)
        {               
            /* Relabel leaves part */
            curr_max += 1;
            for (std::size_t xx=0; xx<sx; ++xx)
                for (std::size_t yy=0; yy<sy; ++yy)
                    if(it3->find(curr_lp[xx][yy])!=it3->end())
                        curr_lp[xx][yy] = curr_max;
            
            /* Copy the labels to the merging sequence */
            curr_nl = curr_nl-it3->size()+1;            
            num_regs_at_th.push_back(pair<double,double>(it->first,curr_nl));
        }
    }
    /*-----------------------------------------------------------*/
    
    
    
    /*-----------------------------------------------------------*/
    /*  Find the threshold at which we have the number of        */
    /*                  regions we want                          */
    /*-----------------------------------------------------------*/
    double curr_nr = 1;
    list<pair<double,double> >::reverse_iterator rit = num_regs_at_th.rbegin();
    if (nreg_max==1)
        ucm_th[0][0] = rit->first;
    else
    {
        double nl_prev = rit->second;
        ++rit;
        for( ; rit!=num_regs_at_th.rend(); ++rit)
        {
            double delta_nl = rit->second-nl_prev;
            curr_nr += delta_nl+1;
            if (curr_nr>nreg_max)
            {
                --rit;
                ucm_th[0][0] = rit->first;
                break;
            }
            nl_prev = rit->second;
        }
    }
    /*-----------------------------------------------------------*/
    
}
