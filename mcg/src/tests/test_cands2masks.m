% ------------------------------------------------------------------------ 
%  Copyright (C)
%  Universitat Politecnica de Catalunya BarcelonaTech (UPC) - Spain
% 
%  Jordi Pont-Tuset <jordi.pont@upc.edu>
%  June 2013
% ------------------------------------------------------------------------ 

%% Dummy dummy
lp = [ 1  1  1
       2  3  2
       2  2  2];
ms_matrix = [ 1  2  4
              3  4  5];

 cands     = [1 2 0
              1 2 3
              2 0 0];
masks = cands2masks(cands, lp, ms_matrix);
assert(isequal(masks(:,:,1),[1 1 1
                             1 0 1
                             1 1 1]))
assert(isequal(masks(:,:,2),[1 1 1
                             1 1 1
                             1 1 1]))
assert(isequal(masks(:,:,3),[0 0 0
                             1 0 1
                             1 1 1]))



%% Real test
im_id = '2008_000008';
load(fullfile(root_dir,'datasets','pascal2012','multiscale','multi', [im_id '.mat']));

curr_hier = ucm2hier(ucm2);
ths{1}.start_ths = curr_hier.start_ths';
ths{1}.end_ths = curr_hier.end_ths';
ms{1} = curr_hier.ms_matrix;
lps = curr_hier.leaves_part;
            
[f_lp,f_ms,cands] = full_cands_from_hiers(lps,ms,ths,[500 500 500 500]');

tic
% Get masks
masks = cands2masks(cands, f_lp, f_ms);
toc

% Compute Areas from candidates
b_feats = compute_base_features(f_lp, f_ms, ucm2);
[~,areas1]  = mex_fast_features(cands-1,b_feats.areas,2);

% Compute areas way 2
areas2 = zeros(size(masks,3),1);
for ii=1:size(masks,3)
    tmp = masks(:,:,ii);
    areas2(ii) = sum(tmp(:));
end

assert(isequal(areas1, areas2))




