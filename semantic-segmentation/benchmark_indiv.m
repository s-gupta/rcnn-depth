function evalRes = benchmark_indiv(outp, imSet, param)
% function evalRes = benchmark_indiv(outp, imSet, param)
%   outp can be a directory containing the images as mats...
%     or it can be a struct with superpixels, and scores
%   outp.imList is the list of images...
%   outp.superpixels
%   outp.scores
%   param.threshPR, param.threshIU, param.thresh, param.imSet

% AUTORIGHTS

  fprintf('Benchmarking SS - %s\n', outp);

  param.imSet = imSet;
  infoFile = param.infoFile;
  ignoreBck = param.ignoreBck;
  fieldName = param.fieldName;

  dt = getMetadata(infoFile);
  nclasses = length(dt.className);
  imList = getImageSet(imSet);
  
  parfor i = 1:length(imList),
    imName = imList{i};
    gtim = getGroundTruth(imName, 'segmentation', infoFile);
    pt = load(fullfile(outp, strcat(imName, '.mat')));
    sc{i} = pt.(fieldName)';
    % if(strcmp(fieldName, 'rawScores')), sc{i} = -sc{i}; end
    gt{i} = accumarray([pt.superpixels(:), gtim(:)+1], 1, [max(pt.superpixels(:)) nclasses+1]);
  end
  
  sc = cat(1, sc{:});
  gt = cat(1, gt{:});
  if(param.ignoreBck)
    gt(:,1) = 0;
  end
  spArea = sum(gt, 2);
  gt = bsxfun(@rdivide, gt(:, 2:end), max(eps, spArea));
  
  for i = 1:nclasses,
    [P(:,i), R(:,i), ap(i,1), ~, ~, ~, iu(:,i)] = calcPR(gt(:,i), sc(:,i), spArea);
  end

  evalRes = struct('P', R, 'R', R, 'ap', ap, 'iu', iu, 'maxIU', max(iu, [], 1)', 'className', {dt.className}, 'param', param);
end
