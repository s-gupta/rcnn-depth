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
function eval_and_save_labels(labels_folder,database,gt_set)

if nargin<2
    database = 'pascal2012';
end
if nargin<3
    gt_set = 'val2012';
end

% Get the name of the folder to refer to the method
if strcmp(labels_folder(end),filesep)
    labels_folder(end) = [];
end
tmp = strfind(labels_folder,filesep);
if isempty(tmp)
    method_name = labels_folder;
    masks_dir = fullfile(root_dir, 'datasets', database, method_name);
else
    method_name = labels_folder(tmp(end)+1:end);
    masks_dir = labels_folder;
end

% Results folder
res_dir   = fullfile(root_dir, 'results', database, method_name);
if ~exist(res_dir,'dir')
    mkdir(res_dir)
end

% Load which images to consider
im_ids = database_ids(database,gt_set);

% Sweep all images
num_images = numel(im_ids);
parfor ii=1:num_images
    curr_id = im_ids{ii};
    res_file = fullfile(res_dir,[curr_id '.mat']);
    
    % Are these candidates already evaluated?
    if ~exist(res_file, 'file')

        % Input file with candidates as labels
        labels_file = fullfile(masks_dir,[curr_id '.mat']);

        % Check if labels are computed
        if ~exist(labels_file, 'file')
            error(['Results ''' labels_file '''not found. Have you computed them?']) 
        end
        
        % Load candidates
        cands = load(labels_file);
        
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
                    
        % Overlap superpixels with GT objects
        n_leaves = length(unique(cands.superpixels));
        superpixels = cands.superpixels;
        superpixels(~valid_pixels) = 0;
        leave_fp   = zeros(n_objs,n_leaves);  % False positives
        leave_int  = zeros(n_objs,n_leaves);  % Intersection with GT
        for kk=1:n_objs
            tmp = hist(double(superpixels(:)).*(gt.object(:)==kk),(0:n_leaves));
            leave_int(kk,:) = tmp(2:end);
            tmp = hist(double(superpixels(:)).*(gt.object(:)~=kk),(0:n_leaves));
            leave_fp(kk,:) = tmp(2:end);
        end
        
        % Create matrix padded with zeros to be compatible with mex_eval_labels
        n_cands    = size(cands.labels,1);
        n_max_labels = length(unique(cands.superpixels));
        label_matrix = zeros(n_cands,n_max_labels);
        for jj=1:n_cands
            label_matrix(jj,1:length(cands.labels{jj})) = cands.labels{jj};
        end
        
        % Compute fp, fn etc. from these values on the superpixels
        inters     = zeros(n_objs,n_cands);
        false_pos  = zeros(n_objs,n_cands);
        false_neg  = zeros(n_objs,n_cands);
        jaccards   = zeros(n_objs,n_cands);
        if n_cands>0
            for kk=1:n_objs
                [inters(kk,:),false_pos(kk,:)] = mex_eval_labels(leave_int(kk,:),leave_fp(kk,:),label_matrix);
                false_neg(kk,:) = true_areas(kk)-inters(kk,:);
                jaccards(kk,:) = inters(kk,:)./(false_pos(kk,:)+true_areas(kk));
            end
        end
        assert(sum(isnan(jaccards(:)))==0)
        
        % Store results
        parsave(res_file,jaccards,inters,false_pos,false_neg,true_areas,obj_classes)
    end
end
end


function parsave(res_file,jaccards,inters,false_pos,false_neg,true_areas,obj_classes) %#ok<INUSD>
    save(res_file, 'jaccards','inters', 'false_pos', 'false_neg','true_areas','obj_classes');
end

