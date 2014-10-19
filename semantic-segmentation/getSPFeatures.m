function [features, superpixels, spArea] = getSPFeatures(imName, paths, param)
  count = 0;
  for j = 1:length(param.featureTyp),
    typ = param.featureTyp{j};
  
    fileName = fullfile(paths.ss_feature_dir, typ, strcat(imName, '.mat'));

    dt = load(fileName);
    selThresh = param.selThresh{j};

    for i = 1:length(selThresh),
      count = count+1;
      F{count} = double(dt.features{selThresh(i)}(:,dt.clusters(:,selThresh(i))));
    end
  end
  features = cat(1,F{:});  
  superpixels = dt.superpixels;
  spArea = accumarray(superpixels(:), 1)';
end
