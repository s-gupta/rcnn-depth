function varargout = edgesEval( model, varargin )
% Run and evaluate structured edge detector on BSDS500.
%
% This function first runs the trained structured edge detector on every
% test or validation image in BSDS then call edgesEvalDir.m to perform the
% actual edge evaluation. edgesEval is specific to the structured edge
% detector and BSDS, edgesEvalDir is general purpose edge evaluation code.
% For example usage of edgesEval see edgesDemo.
%
% USAGE
%  varargout = edgesEval( model, prms )
%
% INPUTS
%  model      - structured edge model trained with edgesTrain
%  prms       - parameters (struct or name/value pairs)
%   .dataType   - ['test'] should be either 'test' or 'val'
%   .name       - [''] name to append to evaluation
%   .opts       - {} list of model opts to overwrite
%   .show       - [0] if true plot results using edgesEvalPlot
%   .pDistr     - [{'type','parfor'}] parameters for fevalDistr
%   .cleanup    - [0] if true delete temporary files
%   .thrs       - [99] number or vector of thresholds for evaluation
%   .maxDist    - [.0075] maximum tolerance for edge match
%
% OUTPUTS
%  varargout  - same outputs as edgesEvalDir
%
% EXAMPLE
%
% See also edgesDemo, edgesDetect, edgesEvalDir, edgesEvalPlot
%
% Structured Edge Detection Toolbox      Version 3.01
% Code written by Piotr Dollar, 2014.
% Licensed under the MSR-LA Full Rights License [see license.txt]

% get default parameters
dfs={'dataType','test', 'name','', 'opts',{}, 'show',0, ...
  'pDistr',{{'type','parfor'}}, 'cleanup',0, 'thrs',99, 'maxDist',.0075, ...
  'test_data', [] };
p=getPrmDflt(varargin,dfs,1);

% load model and update model.opts acoording to opts
if( ischar(model) ), model=load(model); model=model.model; end
for i=1:length(p.opts)/2, model.opts.(p.opts{i*2-1})=p.opts{i*2}; end
rgbd=model.opts.rgbd; model.opts.nms=1;

% get list of relevant directories (image, depth, gt, results)
name = [model.opts.modelFnm,p.name];
%% EDIT BY SG
if(~isempty(p.test_data))
  ids = p.test_data.imlist;
  imgDir = p.test_data.image_dir;
  depDir = p.test_data.dep_dir;
  gtDir = p.test_data.ground_truth_dir;
  ext = p.test_data.ext;
  n = length(ids);
  resDir = fullfile(model.opts.modelDir, p.test_data.dataType, name);
else
  imgDir = fullfile(model.opts.bsdsDir,'images',p.dataType);
  depDir = fullfile(model.opts.bsdsDir,'depth',p.dataType);
  gtDir  = fullfile(model.opts.bsdsDir,'groundTruth',p.dataType);
  ids=dir(imgDir); ids=ids([ids.bytes]>0); ids={ids.name}; n=length(ids);
  ext=ids{1}(end-2:end); for i=1:n, ids{i}=ids{i}(1:end-4); end
  resDir = fullfile(model.opts.modelDir, p.dataType, name);
end

assert(exist(imgDir,'dir')==7); assert(exist(gtDir,'dir')==7);

% run edgesDetect() on every image in imgDir and store in resDir
if( ~exist(fullfile([resDir '-eval'],'eval_bdry.txt'),'file') )
  if(~exist(resDir,'dir')), mkdir(resDir); end
  res=cell(1,n); for i=1:n, res{i}=fullfile(resDir,[ids{i} '.png']); end
  do=false(1,n); for i=1:n, do(i)=~exist(res{i},'file'); end
  ids_do=ids(do); res=res(do); m=length(ids_do);
  parfor i=1:m, 
    id=ids_do{i};
    %
    %% BEGIN EDIT BY SG
    opts = model.opts;
    Iin = imread(fullfile(imgDir,[id '.' ext])); I = Iin;
    D=[]; Din = [];
    ng = zeros(size(I,1), size(I,2), 0);
    if(rgbd), 
      Din=imread(fullfile(depDir,[id '.png']));
      D=single(Din)/1e4; 
    end
    if(rgbd==1), I=D; 
    elseif(rgbd==2), I=cat(3,single(I)/255,D);
    elseif(rgbd==3), 
      I=cat(3,single(I)/255,D,1e-3./D); 
      colorModel = model.colorModel.model;
      cues = getAllCues(Iin, Din, model.opts.rgbd3opts.C, colorModel, ...
        model.opts.rgbd3opts.vars, fullfile_ext(model.opts.rgbd3opts.cueCacheDir,id,'mat')); 
      ng = cat(3, cues{:});
    end
    E=edgesDetect(I,model,ng); imwrite(uint8(E*255),res{i});
    %% END EDIT BY SG
    %
  end
end

% perform actual evaluation using edgesEvalDir
varargout=cell(1,max(1,nargout));
[varargout{:}] = edgesEvalDir('resDir',resDir,'gtDir',gtDir,...
  'pDistr',p.pDistr,'cleanup',p.cleanup,'thrs',p.thrs,'maxDist',p.maxDist,'ids',ids);
if( p.show ), figure(p.show); edgesEvalPlot(resDir,name); end

end
