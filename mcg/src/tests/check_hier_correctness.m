% ------------------------------------------------------------------------ 
%  Copyright (C)
%  Universitat Politecnica de Catalunya BarcelonaTech (UPC) - Spain
% 
%  Jordi Pont-Tuset <jordi.pont@upc.edu>
%  June 2013
% ------------------------------------------------------------------------ 
% 
%  If compute_mcg_cands is run in parallel, since the writing to TXT is not
%  'atomic', sometimes two processes wirte in the same file and
%  thus corrupt it.
%  This script checks the correctness of the hierarchies.
%
% ------------------------------------------------------------------------ 
function check_hier_correctness(params)

if nargin==0
    params = get_params();
end

% Load which images to consider from the params.database (train, val, etc.)
im_ids = database_ids(params.database,params.gt_set_test);

res_dir = fullfile(root_dir,'datasets',params.database,params.mcg_id);

% Sweep all images
num_images = length(im_ids);
disp(['Starting test on: ' res_dir])
for im_id = 1:num_images
    
    hier_file = fullfile(res_dir,[im_ids{im_id} '_hier.txt']);
    
    if ~exist(hier_file,'file')
        disp(['Missing: ' hier_file])
    else
        try
            hier = read_hier(hier_file);
            n_leaves = length(unique(hier.leaves_part));
            n_merges = size(hier.ms_matrix,1);
            if ~isequal(n_merges+n_leaves,hier.ms_matrix(end,end))
                disp(['Bad formed: ' hier_file])
            end
        catch %#ok<CTCH>
            disp(['Bad formed: ' hier_file])
        end
    end
end
disp('Done!')