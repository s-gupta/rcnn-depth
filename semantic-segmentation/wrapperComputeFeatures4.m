function wrapperComputeFeatures4(imSet, ucmThresh, typ, indToDo, detClass, detDir, assignType, pThresh, threshSet, dirNameSuffix, nmsThresh, boxFieldName)
% function wrapperComputeFeatures4(imSet, ucmThresh, typ, indToDo, detClass, detDir, assignType, pThresh, threshSet)
%   imSet - set to compute the features for
%   ucmThesh - what ucm theshold to use for computing these features
%   indToDo - indexes in imSet for which you want to run this thing (too slow)
%   typ - 'self', 'det', 'ssMap', 'detInst', 'detInstV2', 'detInstV4' 
%   detClass - which detectors to use
%   detDir - where is everything from the detectors stored
%   assignTyp - 'maxOverlap' or 'bestScore'
%   pThresh - what precision thresholds to use
%   threshSet - what set ot pick the thresholds on
%   dirNameSuffix - directory in which to store the features
%   nmsThresh - what threshold to use for non-max suppression (uses the exact iu rather than the funny iu from voc-code)
%   boxFieldName - what field in the detector output file is the list of boxes

  % Assume that you have the ucms in the directory ucm, the amodal completions in the directory amodal
  % confMapDir = '/work2/sgupta/rmrc/cachedir/release/cache/outs/svm-full_rmrc_ancFull_tr-train1_val-train2_useVal-1/';
  % detDir = '/work2/sgupta/rmrc/voc-cache/train-all-efficient/singleComp-context-1-extraOctvae-1/';
  % threshSet = 'val1';
  
  % detDir = '/work2/sgupta/rmrc/cachedir/segdetV2/detectors/';
  % threshSet = 'train2';
	
  imList = getImageSet(imSet);
	paths = getPaths();

	spParam.ucmThresh = ucmThresh;

	[x, y] = ndgrid([-1 0 1], [-1 0 1]);
	contextR = [x(:), y(:)]';
  contextParam = struct('sz', 1, 'maxPooling', false, 'contextR', contextR);
  contextParam.fName = @contextFeatures;
  
  switch typ,
    case 'self',
      contextParam.typ = 'self';
      contextParam.confMapDir = confMapDir;
      ds = cell(length(imList), 1);
    case 'ssMap',
      contextParam.typ = 'ssMap';
      contextParam.confMapDir = confMapDir;
      ds = cell(length(imList), 1);
    case 'det',
      %% Load the detections here
      contextParam.typ = 'det';
      cls = detClass;
      for i = 1:length(cls),
        dt = load(fullfile(detDir, sprintf('%s_boxes_%s_2013.mat', cls{i}, imSet)), 'ds');
        dsCat{i} = dt.ds;
        for j = 1:length(imList),
          ds{j}{i} = dsCat{i}{j};
        end
      end
    case 'detInst',
      %% Load the detections here
      contextParam.typ = 'detInst';
      cls = detClass;
      for i = 1:length(cls),
        dt = load(fullfile(detDir, sprintf('%s_boxes_%s_2013.mat', cls{i}, imSet)), 'ds');
        dsCat{i} = dt.ds;
        for j = 1:length(imList),
          ds{j}{i} = dsCat{i}{j};
        end
      end
      contextParam.fName = @detectionFeatures;

    case 'detInstV2',
      %% Load the detections here
      contextParam.typ = 'detInstV2';
      contextParam.assignTyp = assignType; %'maxOverlap'; %bestScore';
      contextParam.pThresh = pThresh;
      cls = detClass;
      for i = 1:length(cls),
        [ex, dsCat{i}, dsThresh(i), recall(i)] = loadDetector(cls{i}, detDir, 'val', imSet, contextParam.pThresh);
        fprintf('.');
        %dt = load(fullfile(detDir, sprintf('%s_boxes_%s_2013.mat', cls{i}, imSet)), 'ds');
        %dsThresh(i) = -1.1;
        %dsCat{i} = dt.ds;
        for j = 1:length(imList),
          ds{j}{i} = dsCat{i}{j};
        end
      end
      contextParam.fName = @detectionFeatures;
      typ = sprintf('detInstV2-%s-%02d', contextParam.assignTyp, round(100*contextParam.pThresh)); 
		  dName = fullfile(paths.featuresDir, typ);
      if(~exist(dName)) mkdir(dName); end

    case 'detInstV4',
      %% Load the detections here
      contextParam.typ = 'detInstV4';
      contextParam.assignTyp = assignType; %'maxOverlap'; %bestScore';
      contextParam.pThresh = pThresh;
      cls = detClass;
      parfor i = 1:length(cls),
        [ex, dsCat{i}, dsThresh(i), recall(i)] = loadDetector(cls{i}, detDir, threshSet, imSet, contextParam.pThresh, nmsThresh, boxFieldName);
      end
      for i = 1:length(cls), for j = 1:length(imList), ds{j}{i} = dsCat{i}{j}; end, end
      
      contextParam.fName = @detectionFeatures;
      typ = sprintf('detInstV4-%s-%02d-%s-nmsThresh%d', contextParam.assignTyp, round(100*contextParam.pThresh), dirNameSuffix, round(100*nmsThresh)); 
		  dName = fullfile(paths.featuresDir, typ);
      if(~exist(dName)) mkdir(dName); end
  end

	%parfor i = 1:length(imList),
  parfor ii = 1:length(indToDo),
    i = indToDo(ii);
		% Load the point cloud
		tt = tic();
		[superpixels, ucm, nSP, spArea] = getSuperpixels(imList{i}, spParam.ucmThresh);
		[clusters] = getAmodalCompletion(imList{i});
  
    Z = getFilledDepthImage(imList{i});

		% Compute generic features
    fprintf('%s - ', imList{i});
		contextData = struct('clusters', clusters, 'superpixels', superpixels, 'imName', imList{i}, 'imsize', size(superpixels), 'ds', {ds{i}}, 'Z', Z, 'thresh', dsThresh);
		
    % Compute Features also saves the features in the cache directory.
		[f, sp2reg] = computeFeatures(imList{i}, paths, typ, contextParam, contextData);
	end
end
