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
function [f_lp, f_ms, cands, start_ths, end_ths] = full_cands_from_hiers(lps,ms,ths,n_cands)

n_hiers  = length(ms);
n_r_cand = size(n_cands,1);
assert(size(n_cands,2)==n_hiers)

% Scan all hierarchies
all_cands = cell(n_hiers,n_r_cand);
for ii=1:n_hiers
    n_r_hier = ms{ii}(end,end);
    assert(length(unique(lps(:,:,ii)))==ms{ii}(1,end)-1)
    assert(ms{ii}(1,end)-1+size(ms{ii},1)==n_r_hier)

    % Get all pairs of neighboring leave regions
    [~, idx_neighbors] = seg2gridbmap(lps(:,:,ii));
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

    % Get the 'n_cands' top candidates from each hierarchy
    % Singletons
    if n_r_hier<=n_cands(1,ii)
        all_cands{ii,1} = (1:n_r_hier)';
    else
        all_cands{ii,1} = (n_r_hier-n_cands(1,ii)+1:n_r_hier)';
    end

%     % Pairs, triplets
%     [all_cands{ii,2}, all_cands{ii,3}] = ...
%         mex_get_tree_cands(double(lps(:,:,ii))-1, double(ms{ii})-1,...
%                              neigh_pairs_min-1, neigh_pairs_max-1,...
%                              [n_cands(2,ii), n_cands(3,ii)]);
                         
    % Pairs, triplets, etc.
    all_cands(ii,2:n_r_cand) = ...
    mex_get_tree_cands(double(lps(:,:,ii))-1, double(ms{ii})-1,...
                         neigh_pairs_min-1, neigh_pairs_max-1,...
                         n_cands(2:end,ii));
                         
end

% Put all them together (for each hierarchy)
full_cands = cell(n_hiers,1);
for ii=1:n_hiers
    % Pre-allocate
    n_tot_cand = 0;
    for jj=1:n_r_cand
        n_tot_cand = n_tot_cand + size(all_cands{ii,jj},1);
    end
    full_cands{ii} = zeros(n_tot_cand,n_r_cand);
    % Store
    curr_n = 0;
    for jj=1:n_r_cand
        full_cands{ii}(curr_n+1:curr_n+size(all_cands{ii,jj},1),1:jj) = all_cands{ii,jj};
        curr_n = curr_n+size(all_cands{ii,jj},1);
    end
end

% Compute the 'unique' hierarchy
lps = double(lps);
luts = cell(n_hiers,1);
for ii=1:n_hiers
    needed_regs = unique(full_cands{ii}(:));
    needed_regs(needed_regs==0) = [];

    [lps(:,:,ii),ms{ii},luts{ii}] = mex_prune_tree_to_regions(lps(:,:,ii)-1,ms{ii}-1,needed_regs-1);
end

% Fuse
[f_lp, f_ms] = fuse_bpts(lps,ms);

% Reduce zeros
f_ms = f_ms(:,logical(sum(f_ms)>0));

% Redo cands (relabel)
n_int = length(unique(f_lp));
curr_n_regs = n_int;
cands = [];
for jj=1:n_hiers
    if ~isempty(ms{jj})
        curr_part = lps(:,:,jj);
        assert(max(curr_part(:))+size(ms{jj},1)==size(luts{jj},1))
        new_cands = zeros(size(full_cands{jj}));
        for xx=1:size(full_cands{jj},1)
            for yy=1:size(full_cands{jj},2)
                if full_cands{jj}(xx,yy)>0
                    pos = logical(full_cands{jj}(xx,yy)==luts{jj}(:,2));
                    if sum(pos)~=1
                        error('Oh oh')
                    end
                    new_cands(xx,yy) = curr_n_regs+luts{jj}(pos,1);
                end
            end
        end
        cands = [cands; new_cands]; %#ok<AGROW>
        curr_n_regs = curr_n_regs + ms{jj}(end,end);
        % Sanity check
        assert(ms{jj}(end,end)==length(unique(lps(:,:,jj))) + size(ms{jj},1))
    end
end
assert(isequal(size(cell2mat(full_cands)), size(cands)))

% Redo ths
start_ths = zeros(1,n_int);
end_ths = zeros(1,n_int);
for jj=1:n_hiers
    if ~isempty(ms{jj})
        curr_part = lps(:,:,jj);
        assert(max(curr_part(:))+size(ms{jj},1)==size(luts{jj},1))

        start_ths = [start_ths ths{jj}.start_ths(luts{jj}(:,2))]; %#ok<AGROW>
        end_ths   = [end_ths ths{jj}.end_ths(luts{jj}(:,2))]; %#ok<AGROW>                
    end
end
assert(length(start_ths)==n_int+size(f_ms,1))


end

