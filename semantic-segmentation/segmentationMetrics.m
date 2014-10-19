function [tp, fp, tn, fn, rawcounts] = segmentationMetrics(prob, gt, threshPR, threshIU) 
  % prob in HxWxnumClass

  numClass = size(prob, 3);

  % Calculate category wise precision and recall
  if(length(threshPR) == 0)
    tp = zeros(0, numClass);
    fp = tp; tn = tp; fn = tp;
  else
    for i = 1:numClass,
      gtI = double(gt == i);
      gtI(gt == 0) = NaN;
      probI = prob(:,:,i);
      [tp(:,i), fp(:,i), tn(:,i), fn(:,i)] = getPR(gtI(:), probI(:), threshPR);
    end
  end

  rawcounts = getRawCounts(gt, prob, threshIU);
end

function rawcounts = getRawCounts(gt, prob, thresh)
  % Calculate Intersection over Union for different thresholds of projection...
  num = size(prob, 3) + 1;
  rawcounts = zeros([num, num, length(thresh)]);
  [conf, ind] = max(prob, [], 3);
  locs = gt >= 0;
  for i = 1:length(thresh),
    resim = ind .* (conf >= thresh(i));
    sumim = 1+gt+resim*num;
    hs = histc(sumim(locs), 1:num*num);
      rawcounts(:,:,i) = rawcounts(:,:,i) + reshape(hs(:), size(rawcounts(:,:,i)));
  end
end

function a = reverse(a)
  a = a(end:-1:1);
end

function [tp fp tn fn] = getPR(gt, out, thresh)
  assert(issorted(thresh));
  out = out+1;
  thresh = thresh+1;
  tp = histc(out(:).*gt(:), thresh);
  tp = reverse(cumsum(reverse(tp)));
  fp = histc(out(:).*(1-gt(:)), thresh);
  fp = reverse(cumsum(reverse(fp)));
  fn = nnz(gt(:) == 1) - tp;
  tn = nnz(gt(:) == 0) - fp;
end
