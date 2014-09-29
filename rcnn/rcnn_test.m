function res = rcnn_test(rcnn_model, imdb, suffix)
% res = rcnn_test(rcnn_model, imdb, suffix)
%   Compute test results using the trained rcnn_model on the
%   image database specified by imdb. Results are saved
%   with an optional suffix.

% AUTORIGHTS
% ---------------------------------------------------------
% Copyright (c) 2014, Ross Girshick
% 
% This file is part of the R-CNN code and is available 
% under the terms of the Simplified BSD License provided in 
% LICENSE. Please retain this notice and LICENSE if you use 
% this file (or any portion of it) in your project.
% ---------------------------------------------------------

conf = rcnn_config();
image_ids = imdb.image_ids;

% assume they are all the same
feat_opts = rcnn_model.feat_opts;
training_opts = rcnn_model.training_opts;
num_classes = length(rcnn_model.classes);

exists_or_mkdir(fullfile(conf.cache_dir, 'detections'));
exists_or_mkdir(fullfile(conf.cache_dir, 'pr-curves'));

try
  aboxes = cell(num_classes, 1);
  for i = 1:num_classes
    load(fullfile_ext(conf.cache_dir, 'detections', [rcnn_model.classes{i} '_boxes_' imdb.name suffix], 'mat'));
    aboxes{i} = boxes;
  end
catch
  aboxes = cell(num_classes, 1);
  box_inds = cell(num_classes, 1);
  for i = 1:num_classes
    aboxes{i} = cell(length(image_ids), 1);
    box_inds{i} = cell(length(image_ids), 1);
  end

  if ~isfield(rcnn_model, 'folds')
    folds{1} = 1:length(image_ids);
  else
    folds = rcnn_model.folds;
  end
  
  fold_for_image = zeros(1, length(imdb.image_ids));
  count = 0;
  for f = 1:length(folds)
    for i = folds{f},
      assert(fold_for_image(i) == 0);
      fold_for_image(i) = f;
    end
  end

  parfor i = 1:length(imdb.image_ids),
    tic_toc_print('%s: test (%s) %d/%d\n', procid(), imdb.name, i, length(image_ids));
    d = get_features(imdb, i, feat_opts, false);
    d.feat = rcnn_scale_features(d.feat, training_opts.feat_norm_mean);
    zsAll{i} = bsxfun(@plus, d.feat*rcnn_model.detectors(fold_for_image(i)).W, rcnn_model.detectors(fold_for_image(i)).B);
    boxesAll{i} = d.boxes;
    gtAll{i} = d.gt;
  end
    
  for i = 1:length(imdb.image_ids),
    boxes = boxesAll{i};
    zs = zsAll{i};
    gt = gtAll{i};
    for j = 1:num_classes,
      z = zs(:,j);
      I = find(~gt);
      aboxes{j}{i} = cat(2, single(boxes(I,:)), z(I));
      allboxes{j}{i} = cat(2, single(boxes), z);
    end
  end
  
  for i = 1:num_classes
    save_file = fullfile_ext(conf.cache_dir, 'detections', [rcnn_model.classes{i} '_boxes_' imdb.name suffix], 'mat');
    parsave(save_file, 'boxes', aboxes{i}, 'allboxes', allboxes{i}, 'imdb', imdb);
  end
end

% ------------------------------------------------------------------------
% Peform AP evaluation
% ------------------------------------------------------------------------
for model_ind = 1:num_classes
  cls = rcnn_model.classes{model_ind};
  res(model_ind) = imdb.eval_func(cls, aboxes{model_ind}, imdb);
  resI = res(model_ind);

  % Generate and save the precision recall curve
  print(res(model_ind).plotHandle, '-djpeg', '-r0', ...
      fullfile(conf.cache_dir, 'pr-curves', [cls '_pr_' imdb.name suffix '.jpg']));

  save(fullfile_ext(conf.cache_dir, 'pr-curves', [cls '_pr_' imdb.name suffix], 'mat'), ...
    '-STRUCT', 'resI');
end

save(fullfile_ext(conf.cache_dir, 'pr-curves', ['results_' imdb.name suffix], 'mat'), 'res', 'imdb');

fprintf('\n~~~~~~~~~~~~~~~~~~~~\n');
fprintf('Results:\n');
aps = [res(:).ap]';
disp(aps);
disp(mean(aps));
fprintf('~~~~~~~~~~~~~~~~~~~~\n');
