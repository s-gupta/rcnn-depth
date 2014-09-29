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
function b_feats = compute_base_features(lp, ms, all_ucms)

n_leaves = length(unique(lp));
n_merges = size(ms,1);
n_regs = n_leaves+n_merges;

%% Areas
b_feats.areas = zeros(n_regs,1);

% Leave areas
tmp = regionprops(lp,'Area');
b_feats.areas(1:n_leaves) = [tmp.Area]';

% Mergings
for jj=1:n_merges
    for kk=1:size(ms,2)-1
        if ms(jj,kk)==0
            break;
        end
        b_feats.areas(n_leaves+jj) = b_feats.areas(n_leaves+jj) + b_feats.areas(ms(jj,kk));
    end
end


%% Intersections
b_feats.intersections = mex_fast_intersections(double(lp)-1, double(ms)-1,b_feats.areas);
for jj=1:length(b_feats.areas)
    b_feats.intersections(jj,jj) = b_feats.areas(jj);
end


%% Bounding boxes
b_feats.bboxes = zeros(n_regs,4);
            
% Compute leaves bboxes
rps = regionprops(lp,'BoundingBox');
bboxes_leaves = ceil(reshape([rps.BoundingBox],4,n_leaves)');
b_feats.bboxes(1:n_leaves,3:4) = bboxes_leaves(:,3:4) + bboxes_leaves(:,1:2)-1;
b_feats.bboxes(1:n_leaves,1:2) = bboxes_leaves(:,1:2);
            
% Evolve through merging sequence
for jj=1:n_merges
    parent = ms(jj,end);
    b_feats.bboxes(parent,:) = b_feats.bboxes(ms(jj,1),:);

    for kk=1:size(ms,2)-1
        if ms(jj,kk)==0
            break;
        end
        son = ms(jj,kk);

        b_feats.bboxes(parent,1) = min(b_feats.bboxes(parent,1),b_feats.bboxes(son,1));
        b_feats.bboxes(parent,2) = min(b_feats.bboxes(parent,2),b_feats.bboxes(son,2));
        b_feats.bboxes(parent,3) = max(b_feats.bboxes(parent,3),b_feats.bboxes(son,3));
        b_feats.bboxes(parent,4) = max(b_feats.bboxes(parent,4),b_feats.bboxes(son,4));
    end
end
     

%% Perimeters and contour sums
[~, idx_neighbors] = seg2gridbmap(lp,1);
K = max(idx_neighbors.matrix_max(:)) + 1;
neigh_pairs_matrix = idx_neighbors.matrix_min+K*idx_neighbors.matrix_max;
neigh_pairs = unique(neigh_pairs_matrix);
neigh_pairs(neigh_pairs==0) = [];
neigh_pairs_min = mod(neigh_pairs,K);
neigh_pairs_max = (neigh_pairs-neigh_pairs_min)/K;

n_leaves = length(unique(lp));
n_merges = size(ms,1);
n_regs = n_leaves+n_merges;
n_max_sons = size(ms,2)-1;


features = {'perimeters','contour_sums'};
for ii=1:length(features)
    
    feature = features{ii};
    if strcmp(feature,'perimeters')
        feat = ones(size(neigh_pairs_matrix));
    elseif strcmp(feature,'contour_sums')
        ucm2 = max(all_ucms,[],3);
        feat = ucm2;
    end

%     % Bottleneck! Implement in MEX!
%     % Allocate
%     common_perimeters = zeros(n_regs,n_regs);
%     border_perimeters = zeros(n_regs,1);
%     % Fill leaves pairs
%     for jj=1:length(neigh_pairs_min)
%         if neigh_pairs_min(jj)==0 % With border
%             border_perimeters(neigh_pairs_max(jj)) = sum((neigh_pairs_matrix(:)==neigh_pairs(jj)).*feat(:));
%         else
%             common_perimeters(neigh_pairs_min(jj), neigh_pairs_max(jj)) = sum((neigh_pairs_matrix(:)==neigh_pairs(jj)).*feat(:));
%         end
%     end

    % ------- MEX implementation ----------
    [common_perimeters, border_perimeters] = ...
        mex_base_perimeters(n_leaves,neigh_pairs_min,neigh_pairs_max,neigh_pairs,neigh_pairs_matrix,feat);
    
    % Allocate the rest of regions
    common_perimeters(n_regs,n_regs) = 0;
    border_perimeters(n_regs) = 0;
    % -------------------------------------

    % Evolve border_perimeters
    for jj=1:n_merges
        parent = ms(jj,end);
        for kk=1:n_max_sons
            curr_son = ms(jj,kk);
            if curr_son==0
                break;
            end
            border_perimeters(parent) = border_perimeters(parent)+border_perimeters(curr_son);
        end
    end

    % Store descendants
    descendants = cell(n_regs,1);
    for jj=1:n_leaves
        descendants{jj} = jj;
    end


    % Fill the whole pairs through merging sequence
    common_perimeters = common_perimeters+common_perimeters';
    for jj=1:n_merges
        parent = ms(jj,end);
        descendants{parent} = parent;
        for kk=1:n_max_sons
            curr_son = ms(jj,kk);
            if curr_son==0
                break;
            end
            common_perimeters(:,parent) = common_perimeters(:,parent)+common_perimeters(:,curr_son);

            descendants{parent} = [descendants{parent} descendants{curr_son}];
        end
        common_perimeters(descendants{parent},parent) = 0;
    end
 
 
 
    common_perimeters = triu(common_perimeters);
    common_perimeters = common_perimeters+common_perimeters';
    for jj=1:n_merges
        parent = ms(jj,end);
        descendants{parent} = parent;
        for kk=1:n_max_sons
            curr_son = ms(jj,kk);
            if curr_son==0
                break;
            end
            common_perimeters(:,parent) = common_perimeters(:,parent)+common_perimeters(:,curr_son);

            descendants{parent} = [descendants{parent} descendants{curr_son}];
        end
        common_perimeters(descendants{parent},parent) = 0;
    end
    common_perimeters = triu(common_perimeters);


    % Compute perimeters
    tmp = common_perimeters+ common_perimeters';
    perimeters = zeros(n_regs,1);
    for jj=1:n_regs
        perimeters(jj) = border_perimeters(jj) + sum(tmp(1:n_leaves,jj));
    end

    clear tmp;
            
    if strcmp(feature,'perimeters')
        b_feats.perimeters = perimeters;
        b_feats.border_perimeters = border_perimeters;
        b_feats.common_perimeters = common_perimeters;
    elseif strcmp(feature,'contour_sums')
        b_feats.contour_sums = perimeters;
        b_feats.border_contours = border_perimeters;
        b_feats.common_contours = common_perimeters;
    end
end
