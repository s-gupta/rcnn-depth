function wrapperComputeFeatures4(imSet, typ, indToDo, detClass, detDir, assignType, pThresh, threshSet, nmsThresh, boxFieldName, sp_dir, out_dir)
% function wrapperComputeFeatures4(imSet, typ, indToDo, detClass, detDir, assignType, pThresh, threshSet)
% Input:
%   imSet - set to compute the features for
%   indToDo - indexes in imSet for which you want to run this thing (too slow)
%   typ - 'self', 'det', 'ssMap', 'detInst', 'detInstV2', 'detInstV4' 
%   detClass - which detectors to use
%   detDir - where is everything from the detectors stored
%   assignTyp - 'maxOverlap' or 'bestScore'
%   pThresh - what precision thresholds to use
%   threshSet - what set ot pick the thresholds on
%   nmsThresh - what threshold to use for non-max suppression (uses the exact iu rather than the funny iu from voc-code)
%   boxFieldName - what field in the detector output file is the list of boxes
%   sp_dir - directory which has the superpixels that can be used
%   out_dir - where to store the results of the computation

  imlist = getImageSet(imSet);
	paths = get_paths();

	[x, y] = ndgrid([-1 0 1], [-1 0 1]);
	contextR = [x(:), y(:)]';
  
  switch typ,
    case 'detection-box',
      %% Load the detections here
      detectionParam.typ = 'detection-box';
      detectionParam.assignTyp = assignType;
      detectionParam.pThresh = pThresh;
      cls = detClass;
      parfor i = 1:length(cls),
        [ex, dsCat{i}, dsThresh(i), recall(i)] = load_detector(cls{i}, detDir, threshSet, imSet, detectionParam.pThresh, nmsThresh, boxFieldName);
      end
      for i = 1:length(cls), for j = 1:length(imlist), ds{j}{i} = dsCat{i}{j}; end, end
      
      detectionParam.fName = @detection_features;
  end

	%parfor i = 1:length(imlist),
  parfor ii = 1:length(indToDo),
    i = indToDo(ii);
    imname = imlist{i};
		
		tt = tic();
    % Load the point cloud
    fileName = fullfile_ext(sp_dir, imname, 'mat');
    dt = load(fileName, 'clusters', 'superpixels');
    
    D = getImage(imname, 'depth'); D = double(D)./1000; D = D.*100; Z = D;

		% Compute generic features
    fprintf('%s - ', imlist{i});
		detectionData = struct('clusters', dt.clusters, 'superpixels', dt.superpixels, ...
      'ds', {ds{i}}, 'Z', Z, 'thresh', dsThresh);
		
    % Compute Features also saves the features in the cache directory.
    [f, sp2reg] = compute_ss_features(imlist{i}, out_dir, detectionParam, detectionData);
	end
end
