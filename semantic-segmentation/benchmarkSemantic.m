function evalRes = benchmarkSemantic(outp, infoFile, imSet)
%   outp is a directory containing the images as mats, with a class label for each pixel
% infoFile is what classmapping to use for the benchmarking
% imSet is the imset to benchmark!

  fprintf('Benchmarking SS - %s', outp);

  dt = getMetadata(infoFile);
  nclasses = length(dt.className);
  imList = getImageSet(imSet);
  
  rawcounts = zeros([nclasses+1, nclasses+1]);

  parfor i = 1:length(imList)
    % load the output file here
    dt = load(fullfile(outp, strcat(imList{i}, '.mat')));
    resim = dt.segmentation;
    gtim = getGroundTruth(imList{i}, 'segmentation', infoFile);
    rawCountsAll{i} = getRawCounts(gtim, resim, nclasses);
  end
  
  % Calculate for the supplied single threshold
  rawcounts = zeros(size(rawCountsAll{1}));
  for i = 1:length(rawCountsAll),
    rawcounts = rawcounts + rawCountsAll{i};
  end
  
  rawcounts(1,:)= 0; % Deleting the pixels that are unlabelled in the ground truth
  accuracies = diag(rawcounts)./(max(1,sum(rawcounts,2)+sum(rawcounts,1)'-diag(rawcounts)));
  accuracies = accuracies(2:end);
  freq = sum(rawcounts,2);
  freq = freq(2:end);
  avacc = mean(accuracies);
  fwavacc = (freq'*accuracies)./sum(freq);
  pixacc = sum(diag(rawcounts))./sum(rawcounts(:));
  
  conf = 100*bsxfun(@rdivide, rawcounts, sum(rawcounts,2)+(sum(rawcounts,2)==0));

  className = dt.className; 

  evalRes.accuracies = accuracies;
  evalRes.pixacc = pixacc;
  evalRes.avacc = avacc;
  evalRes.fwavacc = fwavacc;
  evalRes.conf = conf;
  evalRes.rawcounts = rawcounts;
  evalRes.rawCountsAll = rawCountsAll;
  evalRes.imList = imList;
  evalRes.infoFile = infoFile;
  evalRes.outDir = outp;
  evalRes.className = className;
  evalRes.info = 'The first row and column in rawcounts and counts corresponds to the pixels that were unlabelled in the data set and in the output.';
  
  % Write this in a file somewhere?
  fileName = sprintf('%s-%s-results.mat', outp, imSet);
  fprintf('\nSaving the benchmarking results in %s.\n', fileName);
  save(fileName, '-STRUCT', 'evalRes');
end

function rawcounts = getRawCounts(gt, resim, numClass)
  % Calculate the counts for prediction and ground truth..
  num = numClass + 1;
  rawcounts = zeros([num, num]);
  locs = gt(:) >= 0;
  sumim = 1+gt+resim*num;
  hs = histc(sumim(locs), 1:num*num);
  rawcounts = reshape(hs(:), size(rawcounts));;
end

