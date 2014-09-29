args = {};
if strcmp(jobName, 'edges_to_ucms')
  % Uses the precomptued edge cues to do edge detection at multiple scales
  % This edges at multiple scales are used to generate a multiscale UCM.
  p = get_paths();
  modelfile = fullfile_ext(p.contours_model_dir, 'forest', 'modelNyuRgbd-3', 'mat');
  K = 40;
  for i = 1:K,
    args{i} = {sprintf('all_%d_%d_end', i, K), modelfile, p.ucm_dir, p.contours_cues_dir, cropCamera(getCameraParam('color'))};
  end
  jobParam = struct('numThreads', 1, 'codeDir', pwd(), 'preamble', '', 'matlabpoolN', 0, 'globalVars', {{}}, 'fHandle', @wrapper_rgbd_to_ucm, 'numOutputs', 0);
  resourceParam = struct('mem', 8, 'hh', 20, 'numJobs', K, 'ppn', 1, 'nodes', 1, 'logDir', '/work4/sgupta/pbsBatchDir/', 'queue', 'psi', 'notif', false, 'username', 'sgupta', 'headNode', 'zen');
  [jobId jobDir] = jobParallel(jobName, resourceParam, jobParam, args);
  while ~collectJob(jobDir), pause(120); end
end

if strcmp(jobName, 'benchmark_multi_ucm')
  %% Testing code
  p = get_paths();
  c = benchmarkPaths();
  jobParam = struct('numThreads', 1, 'codeDir', pwd(), 'preamble', '', 'matlabpoolN', 1, 'globalVars', {{}}, 'fHandle', @empty, 'numOutputs', 1);
  resourceParam = struct('mem', 2, 'hh', 1, 'numJobs', 40, 'ppn', 1, 'nodes', 1, 'logDir', '/work4/sgupta/pbsBatchDir/', 'queue', 'psi', 'notif', false, 'username', 'sgupta', 'headNode', 'zen');
  pDistrEval={'type', 'psi', 'pLaunch', struct('jobParam', jobParam, 'resourceParam', resourceParam)};
  edgesEvalDir('resDir', fullfile(p.ucm_dir, 'multi-png', filesep()), ...
    'gtDir', fullfile(c.benchmarkGtDir, filesep()), ...
    'pDistr', pDistrEval, 'cleanup', 0, 'thrs', 99, 'maxDist', 0.011, ...
    'ids', getImageSet('test'));
end

if strcmp(jobName, 'pareto')
  p = get_paths();
  params = nyud_params('root_cache_dir', p.cache_dir, 'feature_id', 'depth', 'depth_features', true, 'camera_matrix', cropCamera(getCameraParam('color')));  
  num_regions = 20000;
  pareto_choose_point(params, num_regions);
end

if strcmp(jobName, 'cache-mcg')
  % Cache the base features for MCG
  jobParam = struct('numThreads', 1, 'codeDir', pwd(), 'preamble', '', 'matlabpoolN', 0, 'globalVars', {{}}, 'fHandle', @cache_mcg_features, 'numOutputs', 0);
  resourceParam = struct('mem', 4, 'hh', 20, 'numJobs', 80, 'ppn', 1, 'nodes', 1, 'logDir', '/work4/sgupta/pbsBatchDir/', 'queue', 'psi', 'notif', false, 'username', 'sgupta', 'headNode', 'zen');
  
  params = nyud_params('root_cache_dir', p.cache_dir, 'feature_id', 'depth', 'depth_features', true, 'camera_matrix', cropCamera(getCameraParam('color')));  
  
  imNames = getImageSet('trainval');
  for i = 1:length(imNames), args{i} = {params, [], imNames{i}, params.cacheDir}; end
  
  [jobId jobDir] = jobParallel(jobName, resourceParam, jobParam, args);
  while ~collectJob(jobDir), pause(120); end
end


if strcmp(jobName, 'rank_training'),
  p = get_paths();
  params = nyud_params('root_cache_dir', p.cache_dir, 'feature_id', 'depth', 'depth_features', true, 'camera_matrix', cropCamera(getCameraParam('color')));  
  rank_training(params, params.cacheDir);
end

if strcmp(jobName, 'region-detect')
  p = get_paths();
  params = nyud_params('root_cache_dir', p.cache_dir, 'feature_id', 'depth', 'depth_features', true, 'camera_matrix', cropCamera(getCameraParam('color')));  
  imlist = getImageSet('all');
  rf = loadvar(params.files.trained_classifier,'rf');
  n_cands = loadvar(params.files.pareto_point,'n_cands');
  
  parfor i = 1:length(imlist),
    % Load the UCMs here
    imname = imlist{i};
    ucms = {};
    for ii = 1:length(params.hier_dirs),
      tmp = load(fullfile_ext(params.hier_dirs{ii}, imname, 'mat')); 
      ucms{ii} = tmp.ucm2; 
    end
    D = getImage(imname, 'depth');
    RD = getImage(imname, 'rawdepth');
    mcg_cache_obj = cache_mcg_features(params, ucms, imname, params.cacheDir);
    candidates = compute_mcg_cands(params, rf, n_cands, mcg_cache_obj, D, RD);
    % Save regions here!
    out_file = fullfile_ext(params.files.output_dir, imname, 'mat');
    parsave(out_file, 'scores', candidates.scores, 'bboxes', candidates.bboxes, 'superpixels', candidates.superpixels, 'labels', candidates.labels);
  end
  stats = eval_labels(params.files.output_dir, 'nyud40Obj', 'test');
  save(sprintf('%s-eval.mat', params.files.output_dir), 'stats', 'params'); 
end

if strcmp(jobName, 'show_region_results')
  p = get_paths();
  params = nyud_params('root_cache_dir', p.cache_dir, 'feature_id', 'depth', 'depth_features', true, 'camera_matrix', cropCamera(getCameraParam('color')));  
  method{1} = params.mcg_id; method{2} = 'old-depth'; method{3} = 'old-nodepth';
  % params = nyud_params('feature_id', 'no_depth', 'depth_features', false, 'camera_matrix', cropCamera(getCameraParam('color')));  
  % method{2} = params.mcg_id; 
  legends = {
    ['Our (MCG (RGBD edges,' char(10) 'RGBD feats.)) [RGBD]'], ...
    ['old - depth'], ...
    ['old - nodepth'], ...
  };
  show_region_results(fullfile(p.output_dir, 'regions'), method, legends);
end
