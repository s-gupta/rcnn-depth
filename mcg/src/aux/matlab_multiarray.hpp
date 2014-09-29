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
#ifndef IMAGEPLUS_MATLAB_MULTIARRAY_HPP
#define IMAGEPLUS_MATLAB_MULTIARRAY_HPP

#include "mex.h"
#include "boost/multi_array.hpp"
#include <typeinfo>

template<std::size_t D>
bool dimension_coherence(const mxArray *data)
{
    return mxGetNumberOfDimensions(data)==D;
}

template<typename T>
bool type_coherence(const mxArray *data)
{
    if(mxGetClassID(data)==mxDOUBLE_CLASS)
        return typeid(double)==typeid(T);
    else if(mxGetClassID(data)==mxLOGICAL_CLASS)
        return typeid(bool)==typeid(T);
    else if(mxGetClassID(data)==mxSINGLE_CLASS)
        return typeid(float)==typeid(T);
    else if(mxGetClassID(data)==mxINT8_CLASS)
        return typeid(char)==typeid(T);
    else if(mxGetClassID(data)==mxUINT8_CLASS)
        return typeid(unsigned char)==typeid(T);
    else if(mxGetClassID(data)==mxINT16_CLASS)
        return typeid(short int)==typeid(T);
    else if(mxGetClassID(data)==mxUINT16_CLASS)
        return typeid(unsigned short int)==typeid(T);
    else if(mxGetClassID(data)==mxINT32_CLASS)
        return typeid(int)==typeid(T);
    else if(mxGetClassID(data)==mxUINT32_CLASS)
        return typeid(unsigned int)==typeid(T);
    else if(mxGetClassID(data)==mxINT64_CLASS)
        return typeid(long long)==typeid(T);
    else if(mxGetClassID(data)==mxUINT64_CLASS)
        return typeid(unsigned long long)==typeid(T);
    else
        return false;
}

template<typename T>
class ConstMatlabMultiArray : public boost::const_multi_array_ref<T,2>{
public:
    ConstMatlabMultiArray(const mxArray *data) : boost::const_multi_array_ref<T,2>((T*)mxGetData(data), boost::extents[mxGetM(data)][mxGetN(data)], boost::fortran_storage_order())
    {
        if(!type_coherence<T>(data))
            mexErrMsgTxt("Data types not coherent" );
        if(!dimension_coherence<2>(data))
            mexErrMsgTxt("Data dimensions not coherent");
    }
};

template<typename T>
class MatlabMultiArray : public boost::multi_array_ref<T,2>
{
public:
    MatlabMultiArray(mxArray *data) : boost::multi_array_ref<T,2>((T*)mxGetData(data), boost::extents[mxGetM(data)][mxGetN(data)], boost::fortran_storage_order())
    {    
        if(!type_coherence<T>(data))
            mexErrMsgTxt("Data types not coherent" );
        if(!dimension_coherence<2>(data))
            mexErrMsgTxt("Data dimensions not coherent" );
    }
};

template<typename T>
class ConstMatlabMultiArray3 : public boost::const_multi_array_ref<T,3>{
public:
    ConstMatlabMultiArray3(const mxArray *data) : boost::const_multi_array_ref<T,3>((T*)mxGetData(data), boost::extents[mxGetDimensions(data)[0]][mxGetDimensions(data)[1]][mxGetDimensions(data)[2]], boost::fortran_storage_order())
    {
        if(!type_coherence<T>(data))
            mexErrMsgTxt("Data types not coherent" );
        if(!dimension_coherence<3>(data))
            mexErrMsgTxt("Data dimensions not coherent");
    }
};

template<typename T>
class MatlabMultiArray3 : public boost::multi_array_ref<T,3>
{
public:
    MatlabMultiArray3(mxArray *data) : boost::multi_array_ref<T,3>((T*)mxGetData(data), boost::extents[mxGetDimensions(data)[0]][mxGetDimensions(data)[1]][mxGetDimensions(data)[2]], boost::fortran_storage_order())
    {    
        if(!type_coherence<T>(data))
            mexErrMsgTxt("Data types not coherent" );
        if(!dimension_coherence<3>(data))
            mexErrMsgTxt("Data dimensions not coherent" );
    }
};

#endif
