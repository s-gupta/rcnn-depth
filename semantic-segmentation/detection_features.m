function features = detectionFeatures(superpixels, sp2reg, param, data)
% function features = detectionFeatures(superpixels, sp2reg, param, data)
  [nSP nReg] = size(sp2reg);
  spArea = histc(superpixels(:),1:nSP)';
  regArea = spArea*sp2reg;
  
  switch param.typ,
    case {'detection-box'},
      oC = @(i) detFeature(data.Z, data.ds{i}, data.thresh(i), sp2reg, superpixels, param.assignTyp);
      oContextF = arrayfun(oC, 1:length(data.thresh), 'UniformOutput', false);
      F{1} = cat(1, oContextF{:});
  end

  features = cat(1,F{:});
  assert(length(find(isnan(features))) == 0, 'Nan Features!!');
end

function f = detFeature(Z, ds, thresh, sp2reg, superpixels, assignTyp)
  % For a particular detector, find overlap of all superpixel bounding boxes with the detector outputs
  bb = regionprops(superpixels, 'BoundingBox');
  bb = cat(1,bb(:).BoundingBox);
  bb(:,[3 4]) = bb(:,[1 2])+bb(:,[3 4]);
  for i = 1:size(sp2reg, 2),
    bbR(i,[1 2]) = min(bb(sp2reg(:,i), [1 2]), [], 1);
    bbR(i,[3 4]) = max(bb(sp2reg(:,i), [3 4]), [], 1);
  end

  %% ok so, we have the bounding box for all regions, find the intersection over union with ds > ind
  if(size(ds, 1) > 0)
    ds = ds(ds(:,5) > thresh, :);
  else
    ds = zeros(0,5);
  end
  ds(end+1,:) = [-1 -1 -1 -1 100];
  [iu inter a1 a2] = bboxOverlap(bbR, ds);
  ds(end, :) = [1 1 1 1 100];

  % Use a function to pick the required detection
  switch assignTyp
    case 'maxOverlap',
      [overlap, id] = max(iu, [], 2);
      overlap = overlap(:);
    case 'bestScore',
      scIU = bsxfun(@times, ds(:,5)'+2, double(iu > 0));
      [~, id] = max(scIU, [], 2);
      overlap = iu(sub2ind(size(iu), 1:size(iu,1), id'));
      overlap = overlap(:);
  end
  % For each detectior window ad region bounding box, find the mean and the median depth value
  zzR = getMeanZ(Z, bbR);
  zzD = getMeanZ(Z, ds);

  
  % Compute features for each superpixel based on what you got here
  f1 = ds(id, 5)'+2;
  f1 = max(f1, 0);
  f2 = (linIt(inter(sub2ind(size(inter), 1:size(inter,1), id')))./linIt(a1))';
  f3 = (linIt(inter(sub2ind(size(inter), 1:size(inter,1), id')))./linIt(a2(id)))';
  f4 = double(overlap == 0)';
  f1(f4 == 1) = 0; f2(f4 == 1) = 0; f3(f4 == 1) = 0;
  f5 = log(max(0.1, f1+0.1));
  f6 = log(max(0.1, f2+0.1));
  f7 = log(max(0.1, f3+0.1));
  f8 = f1.*f2;
  f9 = f1.*f3;
  f10 = (f1).^2.*f2;
  f11 = (f1).^2.*f3;
  f13 = zzR;
  f13(:, f4 == 1) = 0;
  f14 = zzD(:, id);
  f14(:, f4 == 1) = 0;
  f = cat(1, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f13, f14);
end

function f = getMeanZ(z, bbox)
  bbox = round(bbox);
  bbox(:, [1 2]) = max(bbox(:, [1 2]), 1);
  bbox(:, [3]) = min(bbox(:, [3]), size(z,2));
  bbox(:, [4]) = min(bbox(:, [4]), size(z,1));
  f = zeros(2, size(bbox,1));
  for i = 1:size(bbox, 1),
    zi = linIt(z(bbox(i,2):bbox(i,4), bbox(i,1):bbox(i,3)));
    f(:,i) = [median(zi), mean(zi)];
  end
end

function [out ind] = pasteBoxes(sz, ds, op)
  out = zeros(sz);
  out(:) = -inf;
  ind = zeros(sz);
  if(numel(ds) == 0) ds = zeros(0,5); end
  ds(:,1:4) = round(ds(:,1:4));
  
  ds(:,1) = max(ds(:,1), 1);
  ds(:,2) = max(ds(:,2), 1);
  ds(:,3) = min(ds(:,3), sz(2));
  ds(:,4) = min(ds(:,4), sz(1));

  h = ds(:,4)-ds(:,2)+1;
  w = ds(:,3)-ds(:,1)+1;
  for i = 1:size(ds,1),
    % mask = ones(h(i), w(i))*ds(i,5);
    % im = putsubarray(sz, mask, ds(i,2), ds(i,4), ds(i,1), ds(i,3), -inf);
    % %im = putsubarray(sz, mask, ds(i,2), ds(i,4), ds(i,1), ds(i,3), -inf);
    % %out = max(out, im);
    subOut = out(ds(i,2):ds(i,4), ds(i,1):ds(i,3));
    subInd = ind(ds(i,2):ds(i,4), ds(i,1):ds(i,3));
    subInd(subOut < ds(i,5)) = i;
    subOut(subOut < ds(i,5)) = ds(i,5);

    out(ds(i,2):ds(i,4), ds(i,1):ds(i,3)) = subOut;
    ind(ds(i,2):ds(i,4), ds(i,1):ds(i,3)) = subInd;
  end
  out(isinf(out)) = -2;
end
