function evalRes = benchmarkSSInternal(outp, imSet, param)
%   outp can be a directory containing the images as mats...
%   or it can be a struct with superpixels, and scores
% outp.imList is the list of images...
% outp.superpixels
% outp.scores
% param.threshPR, param.threshIU, param.thresh, param.imSet
  fprintf('Benchmarking SS - %s', outp);

  param.imSet = imSet;
  infoFile = param.infoFile;
  thresh = param.thresh;
  threshPR = param.threshPR;
  threshIU = param.threshIU;
  ignoreBck = param.ignoreBck;

  dt = getMetadata(infoFile);
  nclasses = length(dt.className);
  imList = getImageSet(imSet);
  
  tp = zeros(length(threshPR), nclasses);
  fp = tp; fn = tp; tn = tp;
  rawcounts = zeros([nclasses+1, nclasses+1, length(threshIU)]);

  parfor i = 1:length(imList)
    prob = getProbCube(outp, imList{i}, nclasses);
    gtim = getGroundTruth(imList{i}, 'segmentation', infoFile);
    [TP(:,:,i), FP(:,:,i), TN(:,:,i), FN(:,:,i), RAWCOUNTS(:,:,:,i)] = segmentationMetrics(prob, gtim, threshPR, threshIU);
    fprintf('.');
    [~, ~, ~, ~, rawCountsAll{i}] = segmentationMetrics(prob, gtim, 0, thresh);
  end
  tp = sum(TP,3);
  fp = sum(FP,3);
  tn = sum(TN,3);
  fn = sum(FN,3);
  rawcounts = sum(RAWCOUNTS, 4);
  fprintf('\n');
  
  %Calculate the classwise APs.
  detailsAP = [];
  if(length(threshPR) > 0)
    P = max(tp,1)./max(tp+fp,1);
    R = max(tp,1)./max(tp+fn,1);
    for i = 1:nclasses,
      ap(i) = VOCap(R(end:-1:1,i), P(end:-1:1,i));
    end
    detailsAP = struct('ap', ap, 'P', P, 'R', R, 'threshPR', threshPR);
  end

  %Calculate the IU over all thresholds..
  detailsIU = [];
  if(length(threshIU) > 0)
    if(ignoreBck)
      rawcounts(1,:,:) = 0;
    end
    for i = 1:length(threshIU)
      rc = rawcounts(:,:,i);
      accuracies(:,i) = diag(rc)./max(sum(rc,2)+sum(rc,1)'-diag(rc),1);
      freq = sum(rc,2);
    end
    accuracies = accuracies(2:end, :);
    freq = freq(2:end, :);
    avacc = mean(accuracies);
    fwavacc = (freq'*accuracies)./sum(freq);
    detailsIU = struct('accuracies', accuracies, 'avacc', avacc, 'fwavacc', fwavacc, 'threshIU', threshIU, 'freq', freq);
  end
  
  % Calculate for the supplied single threshold
  rawcounts = zeros(size(rawCountsAll{1}));
  for i = 1:length(rawCountsAll),
    rawcounts = rawcounts + full(rawCountsAll{i});
  end

  if(ignoreBck)
    rawcounts(1,:) = 0;
  end
  accuracies = diag(rawcounts)./(max(1,sum(rawcounts,2)+sum(rawcounts,1)'-diag(rawcounts)));
  accuracies = accuracies(2:end);
  freq = sum(rawcounts,2);
  freq = freq(2:end);
  avacc = mean(accuracies);
  fwavacc = (freq'*accuracies)./sum(freq);
  
  conf = 100*bsxfun(@rdivide, rawcounts, sum(rawcounts,2)+(sum(rawcounts,2)==0));

  className = dt.className; 

  evalRes.accuracies = accuracies;
  evalRes.avacc = avacc;
  evalRes.fwavacc = fwavacc;
  evalRes.conf = conf;
  evalRes.rawcounts = rawcounts;
  evalRes.rawCountsAll = rawCountsAll;
  evalRes.imList = imList;
  evalRes.outDir = outp;
  evalRes.className = className;
  evalRes.detailsAP = detailsAP;
  evalRes.detailsIU = detailsIU;
  evalRes.param = param;
  evalRes.thresh = param.thresh;
end
