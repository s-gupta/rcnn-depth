#include <math.h>
#include "mex.h"


/*
    paretomember returns the logical Pareto membership of a set of points.

    synopsis:  front = paretofront(objMat)

    created by Yi Cao
    
    y.cao@cranfield.ac.uk
    
    for compiling type 

    mex paretofront.c
    
*/


void paretofront(bool * front, double * M, unsigned int row, unsigned int col);

void mexFunction( int nlhs, mxArray *plhs[] , int nrhs, const mxArray *prhs[] )
{
	bool * front;
	double * M;
	unsigned int row, col;
	const int  *dims;
    
	if(nrhs == 0 || nlhs > 1)
	{
	    printf("\nsynopsis:   front = paretofront(X)");
	    plhs[0]    = mxCreateDoubleMatrix(0 , 0 ,  mxREAL);
	    return;
	}
	
	M = mxGetPr(prhs[0]);
	dims = mxGetDimensions(prhs[0]);
	row = dims[0];
	col = dims[1];
	
	
	
	/* ----- output ----- */

	plhs[0]    = mxCreateLogicalMatrix (row , 1);
	front = (bool *) mxGetPr(plhs[0]);
	
	
	/* main call */
	paretofront(front,  M, row, col);
}

void paretofront(bool * front, double * M, unsigned int row, unsigned int col)
{
    unsigned int t,s,i,j,j1,j2;
    bool *checklist, coldominatedflag;
    
    checklist = (bool *)mxMalloc(row*sizeof(bool));
    for(t = 0; t<row; t++) checklist[t] = true;
    for(s = 0; s<row; s++) {
        t=s;
        if (!checklist[t]) continue;
        checklist[t]=false;
        coldominatedflag=true;
        for(i=t+1;i<row;i++) {
            if (!checklist[i]) continue;
            checklist[i]=false;
            for (j=0,j1=i,j2=t;j<col;j++,j1+=row,j2+=row) {
                if (M[j1] < M[j2]) {
                    checklist[i]=true;
                    break;
                }
            }
            if (!checklist[i]) continue;
            coldominatedflag=false;
            for (j=0,j1=i,j2=t;j<col;j++,j1+=row,j2+=row) {
                if (M[j1] > M[j2]) {
                    coldominatedflag=true;
                    break;
                }
            }
            if (!coldominatedflag) {     //swap active index continue checking
                front[t]=false;
                checklist[i]=false;
                coldominatedflag=true;
                t=i;
            }
        }
        front[t]=coldominatedflag;
        if (t>s) {
            for (i=s+1; i<t; i++) {
                if (!checklist[i]) continue;
                checklist[i]=false;
                for (j=0,j1=i,j2=t;j<col;j++,j1+=row,j2+=row) {
                    if (M[j1] < M[j2]) {
                        checklist[i]=true;
                        break;
                    }
                }
            }
        }
    }
    mxFree(checklist); 
}