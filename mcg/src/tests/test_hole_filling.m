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

 cands     = [1 2];
[cands_hf] = hole_filling(lp, ms_matrix, cands);
assert(isequal(cands_hf,[1 2 3]))


%% Dummy
lp = [ 1  2  3  4  5
       6  7  8  9 10
      11 12 13 14 15
      16 17 18 19 20];
  
    % 21 21 21 21 24
    % 22 25 25 21 24
    % 22 26 26 23 24
    % 23 23 23 23 24
ms_matrix = [ 1  2  3  4  9 21
              6 11  0  0  0 22
             16 17 18 19 14 23
             10 15  5 20  0 24
              7  8  0  0  0 25
             12 13  0  0  0 26
             22 23 24  0  0 27
             26 27  0  0  0 28
             21 25  0  0  0 29
             28 29  0  0  0 30];

 cands     = [21 22 23];
[cands_hf, cands_comp] = hole_filling(lp, ms_matrix, cands);
assert(isequal(cands_hf  ,[21 22 23 25 26]))
assert(isequal(cands_comp,24))


 cands2     = [21 22 23 24 25];
[cands_hf, cands_comp] = hole_filling(lp, ms_matrix, cands2);
assert(isequal(cands_hf  ,[21 22 23 24 25 26]))
assert(isequal(cands_comp,30))

 cands3    = [1 2 3];
[cands_hf, cands_comp] = hole_filling(lp, ms_matrix, cands3);
assert(isequal(cands_hf  ,[1 2 3]))
assert(isequal(cands_comp,[4 9 25 28]))

%% Real test
im_id = '2008_000009';
load(fullfile(root_dir,'datasets','pascal2012','gPb_mUCM','multi', [im_id '.mat']));

curr_hier = ucm2hier(ucm2);
ths{1}.start_ths = curr_hier.start_ths';
ths{1}.end_ths = curr_hier.end_ths';
ms{1} = curr_hier.ms_matrix;
lps = curr_hier.leaves_part;
            
[f_lp,f_ms,cands] = full_cands_from_hiers(lps,ms,ths,[500 500 500 500]');

tic
[cands_hf, cands_comp] = hole_filling(double(f_lp), double(f_ms), cands);
toc

tic
% Get masks
masks    = cands2masks(cands, f_lp, f_ms);
masks_hf2= false(size(masks));

% Perform hole filling by morphology
for ii=1:size(masks,3)
    masks_hf2(:,:,ii) = (imfill(masks(:,:,ii),'holes')>0);
end
toc

% Is equal?
masks_hf = cands2masks(cands_hf, f_lp, f_ms);
assert(isequal(masks_hf,masks_hf2))

% Was there any hole?
assert(~isequal(masks_hf,masks))

% Check complementaries
n_regs = f_ms(end,end);
masks_comp = cands2masks(cands_comp, f_lp, f_ms);
for ii=1:size(masks,3)
    if isequal(cands_comp(ii),n_regs)
        assert(unique(masks_hf(:,:,ii))==1)
    else
        tmp = double(masks_comp(:,:,ii))+double(masks_hf(:,:,ii));
        assert(unique(tmp)==1)
    end
end

%% Show some results
% id = 2000;
% figure;
% subplot(1,3,1);imshow(masks(:,:,id)>0)
% subplot(1,3,2);imshow(masks_hf(:,:,id)>0)
% subplot(1,3,3);imshow(masks_comp(:,:,id)>0)

