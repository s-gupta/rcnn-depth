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
function [cands_hf, cands_comp] = hole_filling(lp, ms_matrix,cands)

% Get all pairs of neighboring leave regions (including border)
[~, idx_neighbors] = seg2gridbmap(lp,1);
K = max(idx_neighbors.matrix_max(:)) + 1;
neigh_pairs = unique(idx_neighbors.matrix_min+K*idx_neighbors.matrix_max);
neigh_pairs(neigh_pairs==0) = [];
neigh_pairs_min = mod(neigh_pairs,K);
neigh_pairs_max = (neigh_pairs-neigh_pairs_min)/K;

if isrow(neigh_pairs_min)
    neigh_pairs_min = neigh_pairs_min';
end
if isrow(neigh_pairs_max)
    neigh_pairs_max = neigh_pairs_max';
end

[cands_hf, cands_comp] = mex_hole_filling(lp, ms_matrix,cands,neigh_pairs_min,neigh_pairs_max);

