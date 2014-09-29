function d = get_features_i(imdb, image_ids_i, opts, onlyGt)
% function d = get_features_i(imdb, image_ids_i, opts, onlyGt)

% AUTORIGHTS

  d = load(fullfile(opts.featDir, imdb.image_ids{image_ids_i}));

  % % Inject code to sample features from the classes is class_id is not empty
  if(onlyGt)
    sel = find(d.gt);
    assert(isequal(sel, [1:length(sel)]'));
    d.gt = d.gt(sel);
    d.overlap = d.overlap(sel,:);
    d.boxes = d.boxes(sel, :);
    d.feat = d.feat(sel, :);
    d.class = d.class(sel, :);
  end

  d.boxes = double(d.boxes);
  % d.feat = rcnn_pool5_to_fcX(d.feat, opts.param.layer, opts.param.cnnWeightsFile); 
  
end
