function [clusters, Z] = amodal_completionFast(superpixels, pc, thresh)
% function [clusters, Z] = amodal_completionFast(superpixels, pc, thresh)
  %Load ground truth and features

  count = 0;
  spLocal = superpixels;
  nSP = max(superpixels(:));
  
  s  = regionprops(superpixels, 'centroid');
  centroids = round(cat(1, s.Centroid));

  clusters = zeros(nSP, length(thresh));

  iii = 1;
  cache.NR = zeros(3,0);
  cache.b = zeros(1,0);
  cache.nSamples = zeros(2,0);
  cache.distSum = zeros(0,0);
  cache.distCount = zeros(0,0);
  [affinity, cache] = affinityHelper(pc, spLocal, cache, 1:nSP, []);
  Z = zeros(0,3);

  % fprintf('Strating...');
  for i = 1:1000,
    if(iii > length(thresh))
      break;
    end
    [ii jj val] = selectToMerge(affinity, spLocal);
    if(isnan(val)) 
      break;
    end
    if(val > thresh(iii))
      % fprintf('%d, %0.3g;   ', iii, val);
      T = accumarray(superpixels(:),spLocal(:),[],@unique);
      T = make_continuous(T,1);
      clusters(:,iii:end) = repmat(T,1,length(thresh)-iii+1);
      iii = iii+1;
      continue;
    end
    
    count = count+1;
    Z(count,:) = [ii jj val];
    spLocal(spLocal == ii) = count+nSP;
    spLocal(spLocal == jj) = count+nSP;
    [affinity, cache] = affinityHelper(pc, spLocal, cache, count+nSP, [ii jj]);
    % fprintf('.');
  end
  %clusters(:,end+1) = 1;
  %thresh(end+1) = inf;
end

function [ii jj val] = selectToMerge(Affinity, spLocal)
  w = [1/1, 1/2];
  affinity = sum(bsxfun(@times, Affinity(:,:,[4, 1]), reshape(w, [1 1 2])),3);
  affinity(find(eye(size(affinity)))) = NaN;
  spArea = accumarray(spLocal(:),1); 
  affinity(spArea < 1000,:) = NaN;
  affinity(:,spArea < 1000) = NaN;
  [val ind] = min(affinity(:));
  [ii jj] = ind2sub(size(affinity),ind);
  % fprintf('!');
end

function [affinity, cache] = affinityHelper(pc, superpixels, cache, II, toDel)
  %Angle between planes fit to superpixels.
  %Absolute difference between b's for the plane fits.
  %Planarity of the regions
  NR = cache.NR;
  b = cache.b;
  nSamples = cache.nSamples;
  distSum = cache.distSum;
  distCount = cache.distCount;
  
  % Recompute the normals to these things...
  nSP = max(superpixels(:));
  nR = length(II);
  vec = zeros([nSP, nR]);
  vec(sub2ind([nSP, nR], II, 1:nR)) = 1;
  [NR b, ~, ~, ~, nSamples] = compute_region_normals(pc, superpixels, vec);
  NR = [cache.NR, NR];
  b = [cache.b, b];
  nSamples = [cache.nSamples, nSamples];
  nSamples(2, toDel) = 0;

  ind = nSamples(2,:) < 100;
  NR(:,ind) = NaN;
  b(:,ind) = NaN;
  
  affinity{1}(:,:,1) = acosd(abs(1-squareform(pdist(NR','cosine'))));
  affinity{1}(:,:,2) = squareform(pdist(b','euclidean'));
  
  %Jitendra's suggestions
  %Fit a plane to the patch 1 and use it to find the error on the other patch and vice versa.
  for i = II,
    dist = abs(sum(bsxfun(@times,pc,reshape(NR(:,i),[1 1 3])),3) + b(i));
    ind = ~isnan(dist);
    distSum(i,1:nSP) = accumarray(superpixels(ind),dist(ind),[nSP 1],@sum, 0);
    distCount(i,1:nSP) = accumarray(superpixels(ind),~isnan(dist(ind)),[nSP 1],@sum, 0);
  end
  if(numel(toDel) > 0)
    distSum(:,II) = distSum(:,II) + sum(distSum(:, toDel), 2);
    distCount(:,II) = distCount(:,II) + sum(distCount(:, toDel), 2);
  end
  dist1 = distSum./distCount;
  dist2 = dist1';
  affinity{1}(:,:,3) = (dist1 + dist2)./2;
  affinity{1}(:,:,4) = (distSum + distSum')./(distCount+distCount');
  affinity{1}(:,:,5) = min(dist1,dist2);
  affinity{1}(:,:,6) = max(dist1,dist2);

  affinity = cat(3,affinity{:});
  ind = nSamples(2,:) < 100;
  affinity(ind,:,:) = NaN;
  affinity(:,ind,:) = NaN;

  cache.NR = NR;
  cache.b = b;
  cache.nSamples = nSamples;
  cache.distSum = distSum;
  cache.distCount = distCount;
end
