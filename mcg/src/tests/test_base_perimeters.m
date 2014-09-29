% ------------------------------------------------------------------------ 
%  Copyright (C)
%  Universitat Politecnica de Catalunya BarcelonaTech (UPC) - Spain
% 
%  Jordi Pont-Tuset <jordi.pont@upc.edu>
%  June 2013
% ------------------------------------------------------------------------ 

lp = [1 1 2
      1 2 2
      3 4 5];
[~,neighs] = seg2gridbmap(lp,1);
nregs = max(lp(:));
K = max(neighs.matrix_max(:)) + 1;
neigh_pairs_matrix = neighs.matrix_min+K*neighs.matrix_max;
neigh_pairs = unique(neigh_pairs_matrix);
neigh_pairs(neigh_pairs==0) = [];
neigh_pairs_min = mod(neigh_pairs,K);
neigh_pairs_max = (neigh_pairs-neigh_pairs_min)/K;

[perims,b_perims] = mex_base_perimeters(nregs,neigh_pairs_min,neigh_pairs_max,neigh_pairs,neigh_pairs_matrix,ones(size(neigh_pairs_matrix)));
assert(isequal(perims,[0     3     1     0     0
                       0     0     0     1     1
                       0     0     0     1     0
                       0     0     0     0     1
                       0     0     0     0     0]));
                   
assert(isequal(b_perims,[4 3 2 1 2]'));