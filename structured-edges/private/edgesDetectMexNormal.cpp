/*******************************************************************************
* Structured Edge Detection Toolbox      Version 1.00
* Copyright 2013 Piotr Dollar.  [pdollar-at-microsoft.com]
* Please email me if you find bugs, or have suggestions or questions!
* Licensed under the MSR-LA Full Rights License [see license.txt]
*******************************************************************************/
#include <mex.h>
#include <math.h>
#include <stdlib.h>
#ifdef USEOMP
#include <omp.h>
#endif

typedef unsigned int uint32;
typedef unsigned short uint16;
#define min(x,y) ((x) < (y) ? (x) : (y))


// construct lookup array for mapping fids to channel indices
uint32* buildLookup( int *dims, int w ) {
  int c, r, z, n=w*w*dims[2]; uint32 *cids=new uint32[n]; 
  // mexPrintf("Size in build lookup: %d\n", n);
  n=0;
  for(z=0; z<dims[2]; z++) for(c=0; c<w; c++) for(r=0; r<w; r++)
    cids[n++] = z*dims[0]*dims[1] + c*dims[0] + r;
  return cids;
}

// construct lookup arrays for mapping fids for self-similarity channel
void buildLookupSs( uint32 *&cids1, uint32 *&cids2, int *dims, int w, int m ) {
  int i, j, z, z1, c, r; int locs[1024];
  int m2=m*m, n=m2*(m2-1)/2*dims[2], s=int(w/m/2.0+.5);
  cids1 = new uint32[n]; cids2 = new uint32[n]; n=0;
  for(i=0; i<m; i++) locs[i]=uint32((i+1)*(w+2*s-1)/(m+1.0)-s+.5);
  for(z=0; z<dims[2]; z++) for(i=0; i<m2; i++) for(j=i+1; j<m2; j++) {
    z1=z*dims[0]*dims[1]; n++;
    r=i%m; c=(i-r)/m; cids1[n-1]= z1 + locs[c]*dims[0] + locs[r];
    r=j%m; c=(j-r)/m; cids2[n-1]= z1 + locs[c]*dims[0] + locs[r];
  }
}

// [E,ind] = mexFunction(model,chns,chnsSsi, chnsNormal) - helper for edgesDetect.m
void mexFunction( int nl, mxArray *pl[], int nr, const mxArray *pr[] )
{
  // get inputs
  mxArray *model = (mxArray*) pr[0];
  float *chns = (float*) mxGetData(pr[1]);
  float *chnsSs = (float*) mxGetData(pr[2]);
  float *chnsNormal = (float*) mxGetData(pr[3]);

  // extract relevant fields from model and options
  float *thrs = (float*) mxGetData(mxGetField(model,0,"thrs"));
  uint32 *fids = (uint32*) mxGetData(mxGetField(model,0,"fids"));
  uint32 *child = (uint32*) mxGetData(mxGetField(model,0,"child"));
  uint16 *eBins = (uint16*) mxGetData(mxGetField(model,0,"eBins"));
  uint32 *eBnds = (uint32*) mxGetData(mxGetField(model,0,"eBnds"));
  mxArray *opts = mxGetField(model,0,"opts");
  const int shrink = (int) mxGetScalar(mxGetField(opts,0,"shrink"));
  const int imWidth = (int) mxGetScalar(mxGetField(opts,0,"imWidth"));
  const int gtWidth = (int) mxGetScalar(mxGetField(opts,0,"gtWidth"));
  const int nChns = (int) mxGetScalar(mxGetField(opts,0,"nChns"));
  const int nCells = (int) mxGetScalar(mxGetField(opts,0,"nCells"));
  
  const int nNormalCells = (int) mxGetScalar(mxGetField(opts,0,"nNormalCells"));
  const int nNormalFtrs = nNormalCells*nNormalCells*(nNormalCells*nNormalCells-1)/2;
  
  const uint32 nChnFtrs = (uint32) mxGetScalar(mxGetField(opts,0,"nChnFtrs"));
  const uint32 nSimFtrs = (uint32) mxGetScalar(mxGetField(opts,0,"nSimFtrs"));
  const uint32 nNormalFtrsAll = (uint32) mxGetScalar(mxGetField(opts,0,"nNormalFtrs"));
  const int nEdgeBins = (int) mxGetScalar(mxGetField(opts,0,"nEdgeBins"));
  const int stride = (int) mxGetScalar(mxGetField(opts,0,"stride"));
  const int nTreesEval = (int) mxGetScalar(mxGetField(opts,0,"nTreesEval"));
  int nThreads = (int) mxGetScalar(mxGetField(opts,0,"nThreads"));

  // get dimensions and constants
  const mwSize *chnsSize = mxGetDimensions(pr[1]);
  const int h = (int) chnsSize[0]*shrink;
  const int w = (int) chnsSize[1]*shrink;
  const mwSize *fidsSize = mxGetDimensions(mxGetField(model,0,"fids"));
  const int nTreeNodes = (int) fidsSize[0];
  const int nTrees = (int) fidsSize[1];
  const int h1 = (int) ceil(double(h-imWidth)/stride);
  const int w1 = (int) ceil(double(w-imWidth)/stride);
  const int h2 = h1*stride+gtWidth;
  const int w2 = w1*stride+gtWidth;
  const int chnDims[3] = {h/shrink,w/shrink,nChns};
  const int normalDims[3] = {h/shrink,w/shrink,nNormalFtrs};
  const int indDims[3] = {h1,w1,nTreesEval};
  const int outDims[3] = {h2,w2,nEdgeBins};

  // construct lookup tables
  uint32 *eids, *cids, *cids1, *cids2;
  eids = buildLookup( (int*)outDims, gtWidth );
  cids = buildLookup( (int*)chnDims, imWidth/shrink );
  buildLookupSs( cids1, cids2, (int*)chnDims, imWidth/shrink, nCells );
  uint32 *nids;
  nids = buildLookup( (int*)normalDims, 1);

  // mexPrintf("Computed all look up tables.\n");
  // mexPrintf("Size of normalDims - %d, %d, %d\n", normalDims[0], normalDims[1], normalDims[2]);


  // create outputs
  pl[0] = mxCreateNumericArray(3,outDims,mxSINGLE_CLASS,mxREAL);
  float *E = (float*) mxGetData(pl[0]);
  pl[1] = mxCreateNumericArray(3,indDims,mxUINT32_CLASS,mxREAL);
  uint32 *ind = (uint32*) mxGetData(pl[1]);

  // apply forest to all patches and store leaf inds
  // #ifdef USEOMP
  // nThreads = min(nThreads,omp_get_max_threads());
  // #pragma omp parallel for num_threads(nThreads)
  // #endif
  for( int c=0; c<w1; c++ ) for( int t=0; t<nTreesEval; t++ ) {
    for( int r0=0; r0<2; r0++ ) for( int r=r0; r<h1; r+=2 ) {
      int o = (r*stride/shrink) + (c*stride/shrink)*h/shrink;
      // select tree to evaluate
      int t1 = ((r+c)%2*nTreesEval+t)%nTrees; uint32 k = t1*nTreeNodes;
      while( child[k] ) {
        // compute feature (either channel or self-similarity feature)
        uint32 f = fids[k]; float ftr;
        if( f<nChnFtrs ) 
          ftr = chns[cids[f]+o]; 
        else {
          if(f<nChnFtrs+nSimFtrs){
            ftr = chnsSs[cids1[f-nChnFtrs]+o]-chnsSs[cids2[f-nChnFtrs]+o];
          }
          else{
            // for(int i = 0; i < 300; i++) mexPrintf("%d: %d\n", i, nids[i]); mexPrintf("\n\n\n");
            // normal self similarity features
            // mexPrintf("%d, ", f);
            uint32 ff = f-nChnFtrs-nSimFtrs;
            while(ff >= nNormalFtrs) ff = ff-nNormalFtrs;
            // mexPrintf("%d, nids: %d, ", ff, nids[ff]);
            ftr = chnsNormal[nids[ff]+o];
          }
        }
        // compare ftr to threshold and move left or right accordingly
        if( ftr < thrs[k] ) k = child[k]-1; else k = child[k];
        k += t1*nTreeNodes;
      }
      // store leaf index and update edge maps
      ind[ r + c*h1 + t*h1*w1 ] = k;
    }
  }

  // compute edge maps (avoiding collisions from parallel executions)
  for( int c0=0; c0<gtWidth/stride; c0++ ) {
    #ifdef USEOMP
    #pragma omp parallel for num_threads(nThreads)
    #endif
    for( int c=c0; c<w1; c+=gtWidth/stride ) {
      for( int r=0; r<h1; r++ ) for( int t=0; t<nTreesEval; t++ ) {
        uint32 k = ind[ r + c*h1 + t*h1*w1 ];
        float *E1 = E + (r*stride) + (c*stride)*h2;
        int b0=eBnds[k], b1=eBnds[k+1]; if(b0==b1) continue;
        for( int b=b0; b<b1; b++ ) E1[eids[eBins[b]]]++;
      }
    }
  }

  delete [] eids; delete [] cids; delete [] cids1; delete [] cids2; delete [] nids;
}
