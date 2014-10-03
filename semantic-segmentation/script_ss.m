p = get_paths(); imname = 'img_5001';
ucm_thresh = 0.16;
amodal_thresh = [-1 26];
C = cropCamera(getCameraParam('color'));

dt = load(fullfile_ext(p.ucm_dir, 'multi', imname, 'mat'));ucm2 = dt.ucm2;
sp = bwlabel(dt.ucm2 < ucm_thresh); sp = sp(2:2:end, 2:2:end);

RD = getImage(imname, 'rawdepth');
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

% Compute the normals for this image
[N1 b1] = computeNormalsSquareSupport(z./100, missingMask, normalParam.patchSize(1),...
  1, C, ones(size(z)));
[N2 b2] = computeNormalsSquareSupport(z./100, missingMask, normalParam.patchSize(2),...
  1, C, ones(size(z)));
N = N1; 

% Compute the direction of gravity
yDir = getYDir(N2, yDirParam);

genericData = struct('pc', pc, 'clusters', clusters, 'superpixels', sp,...
      'normals', N1, 'yDir', yDir, 'zgMax', zgMax, 'bgOriented', bgOriented, 'bgMax', bgMax);
% Compute Features also saves the features in the cache directory.
[f, sp2reg] = computeFeatures(imList{i}, paths, 'generic', genericParam, genericData);
