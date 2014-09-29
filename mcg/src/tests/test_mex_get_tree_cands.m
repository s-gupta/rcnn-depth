%% Dummy test
lp = [1 1 1 1
      1 2 2 1
      1 3 3 1
      1 1 1 1];
ms = [2 3 4
      1 4 5];

% Get all pairs of neighboring leave regions
[~, idx_neighbors] = seg2gridbmap(lp);
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

n_pairs    = 100;
n_triplets = 100;
curr_cands = mex_get_tree_cands(double(lp)-1, double(ms)-1,...
                                 neigh_pairs_min-1, neigh_pairs_max-1,...
                                 [n_pairs, n_triplets]);
                             
assert(isequal(curr_cands{1},[1 2
                              1 3]))
assert(isempty(curr_cands{2}))                            
                             
%% Simple test 1
lp = [1 2 3 4 5];
ms = [2 3 6
      4 5 7
      6 1 8
      7 8 9];
  
% Get all pairs of neighboring leave regions
[~, idx_neighbors] = seg2gridbmap(lp);
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

n_pairs    = 100;
n_triplets = 100;
curr_cands = mex_get_tree_cands(double(lp)-1, double(ms)-1,...
                                 neigh_pairs_min-1, neigh_pairs_max-1,...
                                 [n_pairs, n_triplets]);
assert(isequal(curr_cands{1},[6 7
                              4 8
                              4 6
                              1 2
                              3 7
                              3 4]))
assert(isempty(curr_cands{2}))

 
%%  Simple test 2
lp = [1 2 3 4];
ms = [1 2 5
      3 4 6
      5 6 7];
  
% Get all pairs of neighboring leave regions
[~, idx_neighbors] = seg2gridbmap(lp);
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

n_pairs    = 100;
n_triplets = 100;
curr_cands = mex_get_tree_cands(double(lp)-1, double(ms)-1,...
                                neigh_pairs_min-1, neigh_pairs_max-1,...
                                [n_pairs, n_triplets]);
assert(isequal(curr_cands{1},[3 5
                              2 6
                              2 3]))
assert(isempty(curr_cands{2}))


%% Simple with 4-tuples
lp = [1 2 3 4 5
      1 3 3 4 5
      1 4 4 4 5
      1 5 5 5 5];
ms = [1 2 6
      3 6 7
      4 7 8
      5 8 9];
  
% Get all pairs of neighboring leave regions
[~, idx_neighbors] = seg2gridbmap(lp);
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

n_cands    = [100, 100, 100, 100];
curr_cands = mex_get_tree_cands(double(lp)-1, double(ms)-1,...
                                neigh_pairs_min-1, neigh_pairs_max-1,...
                                n_cands);
assert(isequal(curr_cands{1},[4 5
                              5 7
                              3 4
                              5 6
                              4 6
                              1 5
                              1 4
                              1 3
                              2 3]))
assert(isequal(curr_cands{2},[3 4 5
                              4 5 6
                              1 4 5
                              1 3 5
                              1 3 4
                              2 3 4]))
assert(isequal(curr_cands{3}, [1 3 4 5
                               2 3 4 5]))
assert(isempty(curr_cands{4}))