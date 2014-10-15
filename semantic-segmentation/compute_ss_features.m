function [features, sp2reg] = compute_ss_features(imname, out_dir, param, data)
% function [features, sp2reg] = compute_ss_features(imname, out_dir, param, data)
  try
    assert(~isempty(imname) && ~isempty(out_dir));
    fileName = fullfile_ext(out_dir, imname, 'mat');
    load(fileName, 'sp2reg', 'features');

  catch
    clusters = data.clusters;
    superpixels = data.superpixels;
    nSP = max(superpixels(:));

    %For each region set compute the regions features 
    for i = 1:size(clusters, 2),
      %calclate the sp2reg for
      cluster = clusters(:,i);
      nR = max(cluster);
      sp2reg{i} = false(nSP,nR);
      sp2reg{i}(sub2ind(size(sp2reg{i}),[1:nSP]',cluster)) = true;
      if(i > 1)
        [TF LOC] = ismember(sp2reg{i}', sp2reg{i-1}', 'rows');
        features{i}(:,find(TF)) = features{i-1}(:,LOC(find(TF)));
        if(~isempty(find(~TF)))
          features{i}(:,find(~TF)) = param.fName(superpixels, sp2reg{i}(:,find(~TF)), param, data);
        end
      else
        %Get features.
        features{i} = param.fName(superpixels, sp2reg{i}, param, data);
      end
      fprintf('.');
    end
    
    features = features;
    sp2reg = sp2reg;

    %Save the things
    if(~isempty(imname) && ~isempty(out_dir))
      fileName = fullfile_ext(out_dir, imname, 'mat');
      save(fileName, 'clusters', 'superpixels', 'sp2reg', 'features');
    end
  end
end
