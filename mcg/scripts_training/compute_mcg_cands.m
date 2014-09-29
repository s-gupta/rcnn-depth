function candidates = compute_mcg_cands(params, rf, n_cands, mcg_cache_obj, D, RD)
% function candidates = compute_mcg_cands(params, rf, n_cands, mcg_cache_obj)

% ------------------------------------------------------------------------ 
%  Copyright (C)
%  Universitat Politecnica de Catalunya BarcelonaTech (UPC) - Spain
%  University of California Berkeley (UCB) - USA
% 
%  Jordi Pont-Tuset <jordi.pont@upc.edu>
%  Pablo Arbelaez <arbelaez@berkeley.edu>
%  June 2014
% ------------------------------------------------------------------------ 
% This file is part of the MCG package presented in:
%    Arbelaez P, Pont-Tuset J, Barron J, Marques F, Malik J,
%    "Multiscale Combinatorial Grouping,"
%    Computer Vision and Pattern Recognition (CVPR) 2014.
% Please consider citing the paper if you use this code.
% ------------------------------------------------------------------------

  n_hiers = length(params.hier_dirs);

  % Load trained random forest
  if(isstr(rf)), rf = loadvar(rf,'rf'); end

  % Load Pareto parameters
  if(isstr(n_cands)), n_cands = loadvar(n_cands,'n_cands'); end

  % Number of regions per candidate
  assert(n_hiers==size(n_cands,2));

  f_lp = mcg_cache_obj.f_lp; f_ms = full(mcg_cache_obj.f_ms); feats = mcg_cache_obj.feats;
  bboxes = mcg_cache_obj.bboxes; red_cands = mcg_cache_obj.red_cands; 
  b_feats_intersections = mcg_cache_obj.b_feats_intersection; 

  if(params.depth_features)
    tt = tic();
    sp = f_lp;
    sp2regC = cands2labels(red_cands, f_ms);
    sp2reg = false(max(sp(:)), length(sp2regC));
    for j = 1:length(sp2regC), sp2reg(sp2regC{j},j) = true; end
    D = double(D)./1000;
    C = params.camera_matrix;
    missingMask = RD == 0;
    fdepth = depthFeatures(sp, sp2reg, D, missingMask, C);
    feats = cat(2, feats, fdepth');
    fprintf('Time for depth features: %0.3f\n', toc(tt));
  end

  % Rank candidates
  class_scores = regRF_predict(feats,rf);
  [scores, ids] = sort(class_scores,'descend');
  red_cands = red_cands(ids,:);
  bboxes = bboxes(ids,:);
  if isrow(scores), scores = scores'; end

  % Max margin
  candidates=[];
  [new_ids, mm_scores] = mex_max_margin(red_cands-1,scores, full(b_feats_intersections), params.theta); %#ok<NASGU>
  cand_labels = red_cands(new_ids,:);
  candidates.scores = scores(new_ids);
  bboxes = bboxes(new_ids,:); 

  % Change the coordinates of bboxes to be coherent with
  % other results from other sources (sel_search, etc.)
  candidates.bboxes = [bboxes(:,2) bboxes(:,1) bboxes(:,4) bboxes(:,3)];

  % Get the labels of leave regions that form each candidates
  candidates.superpixels = f_lp;
  candidates.labels = cands2labels(cand_labels,f_ms);
end
