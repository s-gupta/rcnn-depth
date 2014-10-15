function script_ss(imname)
  p = get_paths('release'); 
  % imname = 'img_5001';
  ucm_thresh = 0.20;
  amodal_thresh = [-1 26];
  C = cropCamera(getCameraParam('color'));

  dt = load(fullfile_ext(p.ucm_dir, 'multi', imname, 'mat'));ucm2 = dt.ucm2;
  sp = bwlabel(dt.ucm2 < ucm_thresh); sp = sp(2:2:end, 2:2:end);

  RD = getImage(imname, 'rawdepth');
  I = getImage(imname, 'images');
  RD = double(RD)./1000; RD(RD == 0) = NaN; RD = RD.*100;
  missingMask = RD == 0;

  [x, y, z] = getPointCloudFromZ(RD, C);
  pc = cat(3, x, y, z);
  tt = tic();
  [clusters, Z] = amodal_completion(sp, pc, amodal_thresh);
  fprintf('Time for amodal completion %0.2f.\n', toc(tt));

  yDirParam.angleThresh = [45 15];
  yDirParam.iter = [5 5];
  yDirParam.y0 = [0 1 0]';
  normalParam.patchSize = [3 10];

  gTextonMapParam.nbins = 30;
  gTextonMapParam.yMinMean = -130;
  gTextonMapParam.yMinOutlier = -90;
  gTextonMapParam.dimensions = 900;
  gTextonMapParam.typ = 'gTexton';

  gTextonBOWParam.dimensions = 900;
  gTextonBOWParam.fName = @bow_features;

  % Load the dictionary
  dt = load('codebook_opp_s1.2_K1000_2000.mat');
  vocab = dt.vocab; 
  siftMapParam.dimensions = size(vocab, 2); 
  siftMapParam.siftParam = struct('siftLib', 'semantic-segmentation/colorDescriptor', 'ds_sampling', 4, 'scales', 1.2, 'descriptor', 'opponentsift');
  siftMapParam.typ = 'colorSift';

  siftBOWParam.fName = @bow_features;
  siftBOWParam.dimensions = siftMapParam.dimensions;

  % Compute the normals for this image
  [N1 b1] = computeNormalsSquareSupport(z./100, missingMask, normalParam.patchSize(1),...
    1, C, ones(size(z)));
  [N2 b2] = computeNormalsSquareSupport(z./100, missingMask, normalParam.patchSize(2),...
    1, C, ones(size(z)));

  % Compute the direction of gravity
  yDir = getYDir(N2, yDirParam);

  genericData = struct('pc', pc, 'clusters', clusters, 'superpixels', sp, 'normals', N1, 'yDir', yDir, ...
    'zgMax', zeros([size(RD), 15]), 'bgOriented', zeros([size(RD), 8]), 'bgMax', zeros([size(RD), 3]));

  genericParam.fName = @generic_features;
  [f, sp2reg] = compute_ss_features(imname, p.ss_generic, genericParam, genericData);

  % Compute sift map and features
  siftData = struct('I', I, 'vocab', vocab);
  [map dimensions] = compute_map(imname, p.ss_map_sift, siftMapParam, siftData);

  siftData = struct('map', map, 'clusters', clusters, 'superpixels', sp);
  [f, sp2reg] = compute_ss_features(imname, p.ss_sift, siftBOWParam, siftData);

  %  % Compute geocentric texton maps and features
  gTextonData = struct('pc', pc, 'N', N1, 'yDir', yDir);
  [map, dimensions] = compute_map([], [], gTextonMapParam, gTextonData);

  gTextonData = struct('map', map, 'clusters', clusters, 'superpixels', sp);
  [f, sp2reg] = compute_ss_features(imname, p.ss_gTexton, gTextonBOWParam, gTextonData);

end
