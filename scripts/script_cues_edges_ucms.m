if strcmp(jobName, 'detect-edge')
  jobParam = struct('numThreads', 1, 'codeDir', pwd(), 'preamble', '', 'matlabpoolN', 1, 'globalVars', {{}}, 'fHandle', @detectEdge, 'numOutputs', 0);
  resourceParam = struct('mem', 2, 'hh', 5, 'numJobs', 40, 'ppn', 1, 'nodes', 1, 'logDir', '/work4/sgupta/pbsBatchDir/', 'queue', 'psi', 'notif', false, 'username', 'sgupta', 'headNode', 'zen')


  inDir = '/work5/sgupta/datasets/nyud2/data/';
  outDir = '/work5/sgupta/eccv14-camera-ready/edges-release/';
  for i = 1:1449,
      args{i} = {fullfile(inDir, 'images', sprintf('img_%04d.png', i+5000)), fullfile(inDir, 'depth', sprintf('img_%04d.png', i+5000)), sprintf('img_%04d', i+5000), 'models-sg/forest/modelNyuRgbd-17-release.mat', [2.0 1.0 0.5], fullfile(outDir, sprintf('img_%04d.mat', i+5000))}; 
  end

  [jobId jobDir] = jobParallel(jobName, resourceParam, jobParam, args);
end


if strcmp(jobName, 'edge-ucm')
  jobParam = struct('numThreads', 1, 'codeDir', pwd(), 'preamble', '', 'matlabpoolN', 1, 'globalVars', {{}}, 'fHandle', @wrapper_rgbd_to_ucm, 'numOutputs', 0);
  resourceParam = struct('mem', 4, 'hh', 5, 'numJobs', 60, 'ppn', 1, 'nodes', 1, 'logDir', '/work4/sgupta/pbsBatchDir/', 'queue', 'psi', 'notif', false, 'username', 'sgupta', 'headNode', 'zen')

  for i = 1:1449,
      args{i} = {sprintf('img_%04d', i+5000)};
  end

  [jobId jobDir] = jobParallel(jobName, resourceParam, jobParam, args);
  while ~collectJob(jobDir), pause(60); end

  jobParam = struct('numThreads', 1, 'codeDir', pwd(), 'preamble', '', 'matlabpoolN', 1, 'globalVars', {{}}, 'fHandle', @empty, 'numOutputs', 1);
  resourceParam = struct('mem', 2, 'hh', 1, 'numJobs', 50, 'ppn', 1, 'nodes', 1, 'logDir', '/work4/sgupta/pbsBatchDir/', 'queue', 'psi', 'notif', false, 'username', 'sgupta', 'headNode', 'psi');
  pDistrEval={'type', 'psi', 'pLaunch', struct('jobParam', jobParam, 'resourceParam', resourceParam)};

  gtDir = '/work4/sgupta/tmp/splitBias/sparse_contour_gradients/nyu_v2/nyuFull/groundTruth/test/';
  resDir = '/work5/sgupta/eccv14-camera-ready/ucms-release/multi-png';
  p.maxDist = 0.011; p.thrs = 99; p.cleanup = 0; p.pDistr = pDistrEval;
  edgesEvalDir('resDir', resDir, 'gtDir', gtDir,...
    'pDistr',p.pDistr,'cleanup',p.cleanup,'thrs',p.thrs,'maxDist',p.maxDist);
end

if strcmp(jobName, 'cache-mcg')
  jobParam = struct('numThreads', 1, 'codeDir', pwd(), 'preamble', '', 'matlabpoolN', 0, 'globalVars', {{}}, 'fHandle', @cache_mcg_features, 'numOutputs', 0);
  resourceParam = struct('mem', 4, 'hh', 20, 'numJobs', 60, 'ppn', 1, 'nodes', 1, 'logDir', '/work4/sgupta/pbsBatchDir/', 'queue', 'psi', 'notif', false, 'username', 'sgupta', 'headNode', 'zen');
  params = nyud_params();
  for i = 1:1449,
      args{i} = {params, sprintf('img_%04d', i+5000)};
  end
  [jobId jobDir] = jobParallel(jobName, resourceParam, jobParam, args);
  while ~collectJob(jobDir), pause(120); end
  rank_training(params);
  compute_mcg_cands(params, getImageSet('train'));
  compute_mcg_cands(params, getImageSet('val'));
  compute_mcg_cands(params, getImageSet('test'));
  stats = eval_labels(labels_folder,database,gt_set);
end

if strcmp(jobName, 'edges_to_ucms')
  % Uses the precomptued edge cues to do multiscale edge detection 
  % which is then used to generate a UCM at multiple scales.
  p = get_paths();
  imlist = getImageSet('all');
  model = load(fullfile_ext(p.contours_model_dir, 'forest', 'modelNyuRgbd-3', 'mat')); 
  model = model.model;
  ucm_dir = p.ucm_dir;
  exists_or_mkdir(fullfile(ucm_dir, 'multi'));
  exists_or_mkdir(fullfile(ucm_dir, 'scale_2.0'));
  exists_or_mkdir(fullfile(ucm_dir, 'scale_1.0'));
  exists_or_mkdir(fullfile(ucm_dir, 'scale_0.5'));
  exists_or_mkdir(fullfile(ucm_dir, 'multi-png'));

  parfor i = 1:length(imlist),
    id = imlist{i};
    I = getImage(id, 'images');
    D = getImage(id, 'depth');
    C = cropCamera(getCameraParam('color'));
    sc = [2 1 0.5];
    cacheFile = fullfile_ext(p.contours_cues_dir, id, 'mat');
    tic; [E, Es, O] = detectEdge(I, D, [], cropCamera(getCameraParam('color')), model, sc, [], cacheFile); toc;
    tic; [ucm2 ucms] = contours_to_ucm(I, sc, Es, O); toc;
    outDir = '/work5/sgupta/eccv14-camera-ready/ucms-release/';
    parsave(fullfile_ext(outDir, 'multi', id, 'mat'), 'ucm2', ucm2);
    parsave(fullfile_ext(outDir, 'scale_2.0', id, 'mat'), 'ucm2', ucms(:,:,1));
    parsave(fullfile_ext(outDir, 'scale_1.0', id, 'mat'), 'ucm2', ucms(:,:,2));
    parsave(fullfile_ext(outDir, 'scale_0.5', id, 'mat'), 'ucm2', ucms(:,:,3));
    imwrite(im2uint8(ucm2(3:2:end, 3:2:end)), fullfile_ext(outDir, 'multi-png', id, 'png'));
  end
end
