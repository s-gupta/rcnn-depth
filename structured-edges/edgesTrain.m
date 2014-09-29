function model = edgesTrain( varargin )
% Train structured edge detector.
%
% For an introductory tutorial please see edgesDemo.m.
%
% USAGE
%  opts = edgesTrain()
%  model = edgesTrain( opts )
%
% INPUTS
%  opts       - parameters (struct or name/value pairs)
%   (1) model parameters:
%   .imWidth    - [32] width of image patches
%   .gtWidth    - [16] width of ground truth patches
%   (2) tree parameters:
%   .nPos       - [5e5] number of positive patches per tree
%   .nNeg       - [5e5] number of negative patches per tree
%   .nImgs      - [inf] maximum number of images to use for training
%   .nTrees     - [8] number of trees in forest to train
%   .fracFtrs   - [1/4] fraction of features to use to train each tree
%   .minCount   - [1] minimum number of data points to allow split
%   .minChild   - [8] minimum number of data points allowed at child nodes
%   .maxDepth   - [64] maximum depth of tree
%   .discretize - ['pca'] options include 'pca' and 'kmeans'
%   .nSamples   - [256] number of samples for clustering structured labels
%   .nClasses   - [2] number of classes (clusters) for binary splits
%   .split      - ['gini'] options include 'gini', 'entropy' and 'twoing'
%   (3) feature parameters:
%   .nOrients   - [4] number of orientations per gradient scale
%   .grdSmooth  - [0] radius for image gradient smoothing (using convTri)
%   .chnSmooth  - [2] radius for reg channel smoothing (using convTri)
%   .simSmooth  - [8] radius for sim channel smoothing (using convTri)
%   .normRad    - [4] gradient normalization radius (see gradientMag)
%   .shrink     - [2] amount to shrink channels
%   .nCells     - [5] number of self similarity cells
%   .rgbd       - [0] 0:RGB, 1:depth, 2:RBG+depth (for NYU data only)
%   (4) detection parameters (can be altered after training):
%   .stride     - [2] stride at which to compute edges
%   .multiscale - [0] if true run multiscale edge detector
%   .sharpen    - [2] sharpening amount (can only decrease after training)
%   .nTreesEval - [4] number of trees to evaluate per location
%   .nThreads   - [4] number of threads for evaluation of trees
%   .nms        - [0] if true apply non-maximum suppression to edges
%   (5) other parameters:
%   .seed       - [1] seed for random stream (for reproducibility)
%   .useParfor  - [0] if true train trees in parallel (memory intensive)
%   .modelDir   - ['models/'] target directory for storing models
%   .modelFnm   - ['model'] model filename
%   .bsdsDir    - ['BSR/BSDS500/data/'] location of BSDS dataset
%
% OUTPUTS
%  model      - trained structured edge detector w the following fields
%   .opts       - input parameters and constants
%   .thrs       - [nNodes x nTrees] threshold corresponding to each fid
%   .fids       - [nNodes x nTrees] feature ids for each node
%   .child      - [nNodes x nTrees] index of child for each node
%   .count      - [nNodes x nTrees] number of data points at each node
%   .depth      - [nNodes x nTrees] depth of each node
%   .eBins      - data structure for storing all node edge maps
%   .eBnds      - data structure for storing all node edge maps
%
% EXAMPLE
%
% See also edgesDemo, edgesChns, edgesDetect, forestTrain
%
% Structured Edge Detection Toolbox      Version 3.01
% Code written by Piotr Dollar, 2014.
% Licensed under the MSR-LA Full Rights License [see license.txt]

% get default parameters
% dfs={'imWidth',32, 'gtWidth',16, 'nPos',5e5, 'nNeg',5e5, 'nImgs',inf, ...
%   'nTrees',8, 'fracFtrs',1/4, 'minCount',1, 'minChild',8, ...
%   'maxDepth',64, 'discretize','pca', 'nSamples',256, 'nClasses',2, ...
%   'split','gini', 'nOrients',4, 'grdSmooth',0, 'chnSmooth',2, ...
%   'simSmooth',8, 'normRad',4, 'shrink',2, 'nCells',5, 'rgbd',0, ...
%   'stride',2, 'multiscale',0, 'sharpen',2, 'nTreesEval',4, ...
%   'nThreads',4, 'nms',0, 'seed',1, 'useParfor',0, 'modelDir','models/', ...
%   'modelFnm','model', 'bsdsDir','BSR/BSDS500/data/', ...
%   'parforOpts', [], ...
%   'ngOpts', struct('ngPath', '/work4/sgupta/tmp/splitBias/sparse_contour_gradients/nyu_v2/nyuFull/ng2/', 'ngVar', 'ng2'), ...
%   };

% Parameters after sweep 3
dfs={'imWidth',16, 'gtWidth',16, 'nEdgeBins',1, 'nPos',8e5, 'nNeg',2e5, ...
  'nImgs',inf, 'nTrees',16, 'fracFtrs',1/4, 'minCount',32, 'minChild',5, ...
  'maxDepth',16, 'discretize','pca', 'nSamples', 16, 'nClasses',8, ...
  'split','gini', 'nOrients',4, 'grdSmooth',1, 'chnSmooth',2, ...
  'simSmooth',4, 'normRad',8, 'shrink',1, 'nCells',5, 'rgbd',3, ...
  'stride',2, 'multiscale',1, 'nTreesEval',4, 'nThreads',4, 'nms',0, ...
  'seed',1, 'useParfor',0, 'modelDir','models/', 'modelFnm','model', ...
  'bsdsDir','BSR/BSDS500/data/', 'trainSet','train', 'parforOpts',[], 'sharpen', 0, ...
  'train_data', [], 'rgbd3opts', struct('C', [], 'vars', {}, 'colorModelFile', '', 'cueCacheDir', '')};

opts = getPrmDflt(varargin,dfs,1);
if(nargin==0), model=opts; return; end

% if forest exists load it and return
% cd(fileparts(mfilename('fullpath')));
forestDir = [opts.modelDir '/forest/'];
forestFn = [forestDir opts.modelFnm];
if(exist([forestFn '.mat'], 'file'))
  load([forestFn '.mat']); return; end

% compute constants and store in opts
nTrees=opts.nTrees; nCells=opts.nCells; shrink=opts.shrink;
opts.nPos=round(opts.nPos); opts.nNeg=round(opts.nNeg);
opts.nTreesEval=min(opts.nTreesEval,nTrees);
opts.stride=max(opts.stride,shrink);
imWidth=opts.imWidth; gtWidth=opts.gtWidth;
imWidth=round(max(gtWidth,imWidth)/shrink/2)*shrink*2;
opts.imWidth=imWidth; opts.gtWidth=gtWidth;
nChnsGrad=(opts.nOrients+1)*2; nChnsColor=3;
if(opts.rgbd==1), nChnsColor=1; end
if(opts.rgbd==2), nChnsGrad=nChnsGrad*2; nChnsColor=nChnsColor+1; end
%%
% BEGIN EDIT BY SG
if(opts.rgbd==3), nChnsGrad=nChnsGrad*3+30+2*(1+5); nChnsColor=nChnsColor+2; end
% END EDIT BY SG
%%
nChns = nChnsGrad+nChnsColor; opts.nChns = nChns;
opts.nChnFtrs = imWidth*imWidth*nChns/shrink/shrink;
opts.nSimFtrs = (nCells*nCells)*(nCells*nCells-1)/2*nChns;
opts.nTotFtrs = opts.nChnFtrs + opts.nSimFtrs; disp(opts);

% generate stream for reproducibility of model
stream=RandStream('mrg32k3a','Seed',opts.seed);

% train nTrees random trees (can be trained with parfor if enough memory)
if(opts.useParfor == 1), 
  parfor i=1:nTrees, trainTree(opts,stream,i); end
elseif(opts.useParfor == 0)
  for i=1:nTrees, trainTree(opts,stream,i); end; 
elseif(opts.useParfor == 2),
  for i=1:nTrees, jobArgs{i} = {opts, stream, i}; end
  tic, s = fevalDistr('trainTree', jobArgs, opts.parforOpts); assert(s == 1); toc;
end

% merge trees and save model
model = mergeTrees( opts );
if(opts.rgbd == 3), model.colorModel = load(opts.rgbd3opts.colorModelFile); end
if(~exist(forestDir,'dir')), mkdir(forestDir); end
save([forestFn '.mat'], 'model', '-v7.3');

end

function model = mergeTrees( opts )
% accumulate trees and merge into final model
nTrees=opts.nTrees; gtWidth=opts.gtWidth;
treeFn = [opts.modelDir '/tree/' opts.modelFnm '_tree'];
for i=1:nTrees
  t=load([treeFn int2str2(i,3) '.mat'],'tree'); t=t.tree;
  if(i==1), trees=t(ones(1,nTrees)); else trees(i)=t; end
end
nNodes=0; for i=1:nTrees, nNodes=max(nNodes,size(trees(i).fids,1)); end
% merge all fields of all trees
model.opts=opts; Z=zeros(nNodes,nTrees,'uint32');
model.thrs=zeros(nNodes,nTrees,'single');
model.fids=Z; model.child=Z; model.count=Z; model.depth=Z;
model.segs=zeros(gtWidth,gtWidth,nNodes,nTrees,'uint8');
for i=1:nTrees, tree=trees(i); nNodes1=size(tree.fids,1);
  model.fids(1:nNodes1,i) = tree.fids;
  model.thrs(1:nNodes1,i) = tree.thrs;
  model.child(1:nNodes1,i) = tree.child;
  model.count(1:nNodes1,i) = tree.count;
  model.depth(1:nNodes1,i) = tree.depth;
  model.segs(:,:,1:nNodes1,i) = tree.hs-1;
end
% remove very small segments (<=5 pixels)
segs=model.segs; nSegs=squeeze(max(max(segs)))+1;
parfor i=1:nTrees*nNodes, m=nSegs(i);
  if(m==1), continue; end; S=segs(:,:,i); del=0;
  for j=1:m, Sj=(S==j-1); if(nnz(Sj)>5), continue; end
    S(Sj)=median(single(S(convTri(single(Sj),1)>0))); del=1; end
  if(del), [~,~,S]=unique(S); S=reshape(S-1,gtWidth,gtWidth);
    segs(:,:,i)=S; nSegs(i)=max(S(:))+1; end
end
model.segs=segs; model.nSegs=nSegs;
% store compact representations of sparse binary edge patches
nBnds=opts.sharpen+1; eBins=cell(nTrees*nNodes,nBnds);
eBnds=zeros(nNodes*nTrees,nBnds);
parfor i=1:nTrees*nNodes
  if(model.child(i) || model.nSegs(i)==1), continue; end %#ok<PFBNS>
  E=gradientMag(single(model.segs(:,:,i)))>.01; E0=0;
  for j=1:nBnds, eBins{i,j}=uint16(find(E & ~E0)'-1); E0=E;
    eBnds(i,j)=length(eBins{i,j}); E=convTri(single(E),1)>.01; end
end
eBins=eBins'; model.eBins=[eBins{:}]';
eBnds=eBnds'; model.eBnds=uint32([0; cumsum(eBnds(:))]);
end
