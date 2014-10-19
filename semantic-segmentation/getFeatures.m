function getFeatures(imName, paths, typ, param)

  switch typ,
    case {'generic', 'categorySpecific', 'scene'},
      % Simply load the files here.
      dt = load(fullfile(paths.featureDir, typ, strcat(imName, '.mat')));

      selThresh = param.selThresh;
      for i = 1:length(selThresh),
        count = count+1;
        F{count} = dt.features{selThresh(i)}(:,dt.clusters(:,selThresh(i)));
      end
      features = cat(1,F{:});  

    case {'colorSift', 'gTexton'},
      % Load the maps and compute bow here, and use the params for computing the histogram for each superpixel
      dt = load(fullfile(paths.featureDir, typ, strcat(imName, '.mat')));
      selThresh = param.selThresh;
      for i = 1:length(selThresh),
        count = count+1;
        F{count} = dt.features{selThresh(i)}(:,dt.clusters(:,selThresh(i)));
      end
      features = cat(1,F{:});  
  end
end


function featureReg = histIt(vq, nbins, sp2reg, superpixels)
  %Quantize the lab into bins from 1 to 10
  featureSP = accumarray([superpixels(vq > 0) vq(vq > 0)], 1, [max(superpixels(:)), nbins])';
  featureReg = featureSP*sp2reg;
  featureReg = sparse(bsxfun(@times, featureReg, 1./sum(featureReg,1)));
end
