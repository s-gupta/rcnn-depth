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
function [leaves_part, merging_sequence] = fuse_bpts(parts,mer_seqs)

[leaves_part,lut] = mex_intersect_hierarchies(uint32(parts));

n_reg_inter = length(lut);
n_hier = size(parts,3);

lut2 = double(cell2mat(lut));
r_ids = (1:n_reg_inter)';

% Get sizes
n_max_sons = 0;
for ii=1:n_hier
    n_max_sons = max(size(mer_seqs{ii},2)-1,n_max_sons);
    n_leaves(ii) = max(lut2(:,ii)); %#ok<AGROW>
    n_merges(ii) = size(mer_seqs{ii},1); %#ok<AGROW>
    n_regs(ii) = n_merges(ii) + n_leaves(ii); %#ok<AGROW>
end

% Get intersection sons forming original leaves
inters_sons = cell(1,n_hier);
for ii=1:n_hier
    inters_sons{ii} = cell(1,n_leaves(ii));
    for jj=1:n_leaves(ii)
        inters_sons{ii}{jj} = r_ids(logical(lut2(:,ii)==jj))';
        n_max_sons = max(n_max_sons, length(inters_sons{ii}{jj}));
    end
end

curr_max_reg = n_reg_inter;

merging_sequence = [];
for ii=1:n_hier
    if ~isempty(mer_seqs{ii})
        % Reproduce original merging sequence
        curr_merges = mer_seqs{ii}(:,1:end-1)+curr_max_reg;
        curr_merges(mer_seqs{ii}(:,1:end-1)==0) = 0;
        new_ms = zeros(n_merges(ii),n_max_sons+1);
        new_ms(:,1:size(curr_merges,2)) = curr_merges;
        new_ms(:,end) = mer_seqs{ii}(:,end)+curr_max_reg;

        % Create artificial merging sequence
        art_ms = zeros(n_leaves(ii),n_max_sons+1);
        for jj=1:n_leaves(ii)
            art_ms(jj,1:length(inters_sons{ii}{jj})) = inters_sons{ii}{jj};
            art_ms(jj,end) = jj+curr_max_reg;
        end

        % Stack everything
        merging_sequence = [merging_sequence; art_ms; new_ms]; %#ok<AGROW>
        curr_max_reg = curr_max_reg+n_regs(ii);
    end
end