if strcmp(jobName, 'compute_edge_cues')
  p = get_paths();
  imlist = getImageSet('all');
  colorModel = load(p.rgb_edge_model);
  colorModel = colorModel.model;
  C = cropCamera(getCameraParam('color'));
  vars = {'ng1', 'ng2', 'Es', 'y1', 'y2', 'y3', 'angl1', 'angl2'};
  cueCacheDir = p.contours_cues_dir;

  for i = 1:length(imlist),
    cacheFile = fullfile_ext(cueCacheDir, imlist{i}, 'mat');
    getAllCues(getImage(imlist{i}, 'images'), getImage(imlist{i}, 'depth'), C, colorModel, vars, cacheFile);
  end
end
  
if strcmp(jobName, 'train_edge_model')
  %% set opts for training (see edgesTrain.m)
  p = get_paths();
  c = benchmarkPaths();
  opts=edgesTrain();                                      % default options (good settings)
  opts.modelDir=p.contours_model_dir;                     % model will be in models/forest
  % opts.nPos=5e5; opts.nNeg=5e5;                         % decrease to speedup training
  opts.useParfor=0;                                       % parallelize if sufficient memory
  opts.fracFtrs=1/8;                                      % sample fewer ftrs since using depth
  opts.rgbd=3;                                           % specify use of rgb+d images
  opts.modelFnm=sprintf('modelNyuRgbd-%d', opts.rgbd);    % model name

  opts.train_data = struct('imlist', {getImageSet('train')}, 'image_dir', fullfile(c.dataDir, 'images', filesep()), ...
    'dep_dir', fullfile(c.dataDir, 'depth', filesep()), ...
    'ground_truth_dir', fullfile(c.benchmarkGtDir, filesep()), 'ext', 'png'); 

  jobParam = struct('numThreads', 1, 'codeDir', pwd(), 'preamble', '', 'matlabpoolN', 1, 'globalVars', {{}}, 'fHandle', @empty, 'numOutputs', 1);
  resourceParam = struct('mem', 46, 'hh', 4, 'numJobs', 4, 'ppn', 1, 'nodes', 1, ...
    'logDir', '/work4/sgupta/pbsBatchDir/', 'queue', 'psi', 'notif', false, 'username', 'sgupta', 'headNode', 'zen');
  pDistrTrain={'type', 'psi', 'pLaunch', struct('jobParam', jobParam, 'resourceParam', resourceParam)};
  opts.useParfor = 2; opts.parforOpts = pDistrTrain;

  if(opts.rgbd == 3)
    rgbd3opts.vars = {'ng1', 'ng2', 'Es', 'y1', 'y2', 'y3', 'angl1', 'angl2'};
    rgbd3opts.colorModelFile = p.rgb_edge_model;
    rgbd3opts.cueCacheDir = p.contours_cues_dir;
    rgbd3opts.C = cropCamera(getCameraParam('color'));
    opts.rgbd3opts = rgbd3opts;
  end

  %% train edge detector (~30m/15Gb per tree, proportional to nPos/nNeg)
  tic, model=edgesTrain(opts); toc; % will load model if already trained
end

if strcmp(jobName, 'test_edge_model')
  p = get_paths();
  c = benchmarkPaths();
  model = load(fullfile_ext(p.contours_model_dir, 'forest', 'modelNyuRgbd-3', 'mat')); model = model.model;
  test_data = struct('dataType', 'test', 'imlist', {getImageSet('test')}, ...
    'image_dir', fullfile(c.dataDir, 'images', filesep()), ...
    'dep_dir', fullfile(c.dataDir, 'depth', filesep()), ...
    'ground_truth_dir', fullfile(c.benchmarkGtDir, filesep()), 'ext', 'png', ...
    'res_dir', fullfile(p.contours_model_dir, 'test', 'modelNyuRgbd-3')); 
  
  %% Testing code
  jobParam = struct('numThreads', 1, 'codeDir', pwd(), 'preamble', '', 'matlabpoolN', 1, 'globalVars', {{}}, 'fHandle', @empty, 'numOutputs', 1);
  resourceParam = struct('mem', 2, 'hh', 1, 'numJobs', 50, 'ppn', 1, 'nodes', 1, 'logDir', '/work4/sgupta/pbsBatchDir/', 'queue', 'psi', 'notif', false, 'username', 'sgupta', 'headNode', 'psi');
  pDistrEval={'type', 'psi', 'pLaunch', struct('jobParam', jobParam, 'resourceParam', resourceParam)};

  %% set detection parameters (can set after training)
  model.opts.multiscale=1;          % for top accuracy set multiscale=1
  model.opts.nTreesEval=4;          % for top speed set nTreesEval=1
  model.opts.nThreads=4;            % max number threads for evaluation
  model.opts.nms=0;                 % set to true to enable nms
  model.opts.sharpen=0;             % for top speed set sharpen=0
  if(1), edgesEval( model, 'show', 1, ...
    'name','_release', ...
    'maxDist',.011, 'pDistr', pDistrEval, ...
    'test_data', test_data);
  end
end
