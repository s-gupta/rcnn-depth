function features = bow_features(superpixels, sp2reg, param, data)
% function features = bow_features(superpixels, sp2reg, param, data)
% Constructs BOW features on sift/gtexton maps
% Input: 
%   param.nbins - number of bins that there are
%   data.map - the map containing perpixel words
%   superpixels and sp2reg define the regions over which to accumulate the histograms

  [nSP nReg] = size(sp2reg);
  spArea = histc(superpixels(:),1:nSP)';
  regArea = spArea*sp2reg;
  
  tic;
  nbins = param.dimensions;
  map = data.map;
  F{1} = histIt(map, nbins, sp2reg, superpixels);
  fprintf('BOW Features : %0.3f seconds.\n',toc);
  
  features = cat(1,F{:});
end


function featureReg = histIt(vq, nbins, sp2reg, superpixels)
  %Quantize the lab into bins from 1 to 10
  featureSP = accumarray([superpixels(vq > 0) vq(vq > 0)], 1, [max(superpixels(:)), nbins])';
  featureReg = featureSP*sp2reg;
  featureReg = sparse(bsxfun(@times, featureReg, 1./sum(featureReg,1)));
end
