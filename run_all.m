function [E, ucm2, candidates, detection_scores_no_nms, cls] = run_all(I, D, RD, C, out_file)
% function [E, ucm2, candidates, detection_scores_no_nms, cls] = run_all(I, D, RD, C, out_file)

% AUTORIGHTS

  %% Compute the UCMs
  p = get_paths();
  edge_model_file = fullfile_ext(p.contours_model_dir, 'forest', 'modelNyuRgbd-3', 'mat');
  model = load(edge_model_file);
  model = model.model;
  sc = [2 1 0.5];
  [E, Es, O] = detectEdge(I, D, [], C, model, sc, [], []);
  [ucm2 ucms] = contours_to_ucm(I, sc, Es, O);
  if(~isempty(out_file)), save(out_file, 'E', 'Es', 'O', 'ucm2', 'ucms'); end

  %% Compute the regions
  params = nyud_params('root_cache_dir', p.cache_dir, 'feature_id', 'depth', 'depth_features', true, 'camera_matrix', C);  
  rf = loadvar(params.files.trained_classifier,'rf');
  n_cands = loadvar(params.files.pareto_point,'n_cands');
   
  mcg_cache_obj = cache_mcg_features(params, {ucm2, ucms(:,:,1), ucms(:,:,2), ucms(:,:,3)}, [], []);
  candidates = compute_mcg_cands(params, rf, n_cands, mcg_cache_obj, D, RD);
  if(~isempty(out_file)), save(out_file, '-append', 'candidates'); end


  % Display the superpixels and the regions
  figure(1); 
  subplot(2,3,1); imagesc(Es{2}); axis image; title('Edge Signal');
  subplot(2,3,2); imagesc(ucm2(3:2:end, 3:2:end)); axis image; title('Multi UCM');
  sp = bwlabel(ucm2 < 0.20); sp = sp(2:2:end, 2:2:end);
  for i = 1:3, csp(:,i) = accumarray(sp(:), linIt(I(:,:,i)), [], @mean); end
  subplot(2,3,3); imagesc(ind2rgb(sp, im2double(uint8(csp)))); axis image; title('Superpixels');
  
  boxes = candidates.bboxes(1:2000, [2 1 4 3]);
  
  % Compute the RGB Features
  net_file = fullfile_ext(p.snapshot_dir, sprintf('nyud2_finetune_color_iter_%d', 30000), 'caffemodel');
  net_def_file = fullfile('nyud2_finetuning', 'imagenet_color_256_fc6.prototxt');
  mean_file = fullfile_ext(p.mean_file_color, 'mat');
  rcnn_model = rcnn_create_model(net_def_file, net_file, mean_file);
  rcnn_model = rcnn_load_model(rcnn_model);
  feat{1} = rcnn_features(I, boxes, rcnn_model);

  % Compute the HHA Features
  HHA = saveHHA([], C, [], D, RD);
  net_file = fullfile_ext(p.snapshot_dir, sprintf('nyud2_finetune_hha_iter_%d', 30000), 'caffemodel');
  net_def_file = fullfile('nyud2_finetuning', 'imagenet_hha_256_fc6.prototxt');
  mean_file = fullfile_ext(p.mean_file_hha, 'mat');
  rcnn_model = rcnn_create_model(net_def_file, net_file, mean_file);
  rcnn_model = rcnn_load_model(rcnn_model);
  feat{2} = rcnn_features(HHA, boxes, rcnn_model);
  
  feat = cat(2, feat{:});
  % Load the detectors!
  global RCNN_CONFIG_OVERRIDE;
  conf_override.sub_dir = sprintf('rgb_hha_%d_%s', 30000, 'trainval');
  RCNN_CONFIG_OVERRIDE = @() conf_override;
  conf = rcnn_config();
  dt = load([conf.cache_dir 'rcnn_model'], 'rcnn_model'); rcnn_model = dt.rcnn_model; clear dt;
  feat = rcnn_scale_features(feat, rcnn_model.training_opts.feat_norm_mean);
  detection_scores_no_nms = bsxfun(@plus, feat*rcnn_model.detectors.W, rcnn_model.detectors.B);
  cls = rcnn_model.classes;
  if(~isempty(out_file)), save(out_file, '-append', 'detection_scores_no_nms', 'cls'); end
  
  % Visualize some detections
  cls_id = [2 5 16 17];
  cols = lines(length(cls_id));
  Idet = I;
  for i = 1:length(cls_id),
    dt = load(fullfile(conf.cache_dir, 'pr-curves', sprintf('%s_pr_nyud2_test_release.mat', cls{cls_id(i)})));
    bbox = cat(2, boxes, detection_scores_no_nms(:,cls_id(i)));
    keep = false(size(bbox(:,1)));
    keep(rcnn_nms(bbox, 0.3)) = 1;
    thresh = dt.thresh(find(dt.prec > 0.8, 1, 'last'));
    ind = bbox(:,5) > thresh;
    keep = find(keep & ind);
    bbox = bbox(keep,:);
    if(size(bbox,1) > 0)
      Idet = draw_rect_vec(Idet, bbox(:,1:4)', im2uint8(cols(i,:)), 2);
    end
  end
  figure(1); subplot(2,3,4); imagesc(Idet); axis image; title(['detections - ', sprintf('%s, ', cls{cls_id})]);
  figure(1); subplot(2,3,5); plot([1:length(cls_id)], 1); legend(cls(cls_id)); axis image;
  
  %% Do instance segmentation
  % [] = instance_segmentation(I, D, detections, sp);
  %% Visualize the instance segmentations
end
