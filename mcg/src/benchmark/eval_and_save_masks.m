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
function eval_and_save_masks(masks_folder,database,gt_set)

if nargin<2
    database = 'pascal2012';
end
if nargin<3
    gt_set = 'val2012';
end

% Get the name of the folder to refer to the method
if strcmp(masks_folder(end),filesep)
    masks_folder(end) = [];
end
tmp = strfind(masks_folder,filesep);
if isempty(tmp)
    method_name = masks_folder;
    masks_dir = fullfile(root_dir, 'datasets', database, method_name);
else
    method_name = masks_folder(tmp(end)+1:end);
    masks_dir = masks_folder;
end

% Results folder
res_dir   = fullfile(root_dir, 'results', database, method_name);
if ~exist(res_dir,'dir')
    mkdir(res_dir)
end

% Load which images to consider
im_ids = database_ids(database,gt_set);

% Sweep all images in parallel
matlabpool(4)
num_images = numel(im_ids);
parfor ii=1:num_images
    curr_id = im_ids{ii};
    res_file = fullfile(res_dir,[curr_id '.mat']);
    
    % Are these masks already evaluated?
    if ~exist(res_file, 'file')

        % Input file, with a variable named 'masks'
        masks_file = fullfile(masks_dir,[curr_id '.mat']);

        % Check if masks are already built
        if ~exist(masks_file, 'file')
            error(['Results ''' masks_file '''not found. Have you computed them?']) 
        end
        
        % Load masks
        masks = loadvar(masks_file,'masks');
        
        % Load GT
        gt = get_ground_truth(database,curr_id);

        % Get all objects ids
        obj_ids = unique(gt.object);
        obj_ids(obj_ids==0) = [];
        obj_ids(obj_ids==255) = [];
        n_objs = numel(obj_ids);
        valid_pixels = (gt.object~=255);

        % Store true_areas
        true_areas = zeros(n_objs,1);
        for kk=1:n_objs
            true_areas(kk) = sum(gt.object(:)==kk);
        end
            
        % Store the class of each object
        obj_classes = zeros(n_objs,1);
        for kk=1:n_objs
            idx = find(gt.object==obj_ids(kk),1,'first');
            obj_classes(kk) = gt.class(idx);
        end
                    
        % Sweep and compare all masks
        [areas,inters,false_neg] = mex_eval_masks(masks,gt.object,valid_pixels);
        false_pos = repmat(areas,n_objs,1)-inters;
        jaccards = inters./(inters+false_pos+false_neg);
        assert(sum(isnan(jaccards(:)))==0)

        % Store results
        parsave(res_file,jaccards,inters,false_pos,false_neg,true_areas,obj_classes)
    end
end

matlabpool close

end


function parsave(res_file,jaccards,inters,false_pos,false_neg,true_areas,obj_classes) %#ok<INUSD>
    save(res_file, 'jaccards','inters', 'false_pos', 'false_neg','true_areas','obj_classes');
end

