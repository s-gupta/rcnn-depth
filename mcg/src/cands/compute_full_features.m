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
function [features, bboxes] = compute_full_features(cands,b_feats)

% full features
[~,bboxes]  = mex_fast_features(cands-1,b_feats.bboxes,4);
curr_end_ths = mex_fast_features(cands-1,b_feats.end_ths,1);
curr_start_ths = mex_fast_features(cands-1,b_feats.start_ths,1);
curr_perims  = mex_fast_features(cands-1,b_feats.perimeters,3,b_feats.common_perimeters);
curr_th_sums = mex_fast_features(cands-1,b_feats.contour_sums,3,b_feats.common_contours);


% Mean contour strength
curr_m_cont = curr_th_sums./curr_perims;

% Areas
[curr_area_bal,curr_areas]  = mex_fast_features(cands-1,b_feats.areas,2);

% Compacity
curr_comp = (curr_perims.*curr_perims)./curr_areas;

% Contour sum / sqrt(area)
curr_cs_area = curr_th_sums./sqrt(curr_areas);

% Bounding box
norm_bboxes = [(bboxes(:,1)-1)/b_feats.im_size(2) (bboxes(:,2)-1)/b_feats.im_size(1) bboxes(:,3)/b_feats.im_size(2) bboxes(:,4)/b_feats.im_size(1)];

area_bboxes = (bboxes(:,3)-bboxes(:,1)+1).*(bboxes(:,4)-bboxes(:,2)+1);
area_per_bbox = curr_areas./area_bboxes;
aspect_ratio = (bboxes(:,3)-bboxes(:,1)+1)./(bboxes(:,4)-bboxes(:,2)+1);

% Append current features
features = [curr_end_ths curr_start_ths curr_area_bal curr_areas/(b_feats.im_size(1)*b_feats.im_size(2)) curr_perims curr_th_sums curr_m_cont curr_comp curr_cs_area norm_bboxes area_bboxes area_per_bbox aspect_ratio];
end

