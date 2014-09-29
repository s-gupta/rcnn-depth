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
function stats = eval_masks(masks_folder,database,gt_set,n_cands)

if nargin<2
    database = 'pascal2012';
end
if nargin<3
    gt_set = 'val2012';
end
if nargin<4
    % Number of sampled candidates
    n_cands = [10:5:100,125:25:1000,1500:500:6000,10000];
end

% Get the name of the folder to refer to the method
if strcmp(masks_folder(end),filesep)
    masks_folder(end) = [];
end
tmp = strfind(masks_folder,filesep);
if isempty(tmp)
    method_name = masks_folder;
else
    method_name = masks_folder(tmp(end)+1:end);
end

res_dir    = fullfile(root_dir,'results',database, method_name);
stats_file = fullfile(root_dir,'results',database, [method_name '_' gt_set '.mat']);

% Is the result already gathered?
if exist(stats_file, 'file')
    load(stats_file) 
    recompute = 0;
    disp(['Loaded: ' stats_file '.'])
else
    disp(['RECOMPUTING: ' stats_file '.'])
    recompute = 1;
end

if recompute
    %% Evaluate and save each image independently to be able to parallelize
    % You can adapt the matlabpool to your system
    eval_and_save_masks(masks_folder,database,gt_set);
    
    %% Gather and save results
    % Load which images to consider
    im_ids = database_ids(database,gt_set);
    
    % Store and initialize
    stats.n_cands = n_cands;
    stats.gt_set  = gt_set;
    stats.obj_classes = [];
    
    % Compute statistics
    stats.num_objects = 0;
    for ii=1:numel(im_ids)
        res_file = fullfile(res_dir, [im_ids{ii} '.mat']);
        if exist(res_file,'file')
            load(res_file)
            for jj=1:length(stats.n_cands)
                curr_n_cands = min(stats.n_cands(jj), size(jaccards,2));
                to_consider = 1:curr_n_cands;
                stats.all_n_masks(ii,jj) = length(to_consider);
                if (stats.all_n_masks(ii,jj)>0) 
                    for kk=1:size(jaccards,1)
                        [stats.max_J(stats.num_objects+kk,jj), which_one] = max(jaccards(kk,to_consider));
                        stats.max_indicator(stats.num_objects+kk,jj) = to_consider(which_one);
                        stats.max_fp(stats.num_objects+kk,jj) = false_pos(kk,to_consider(which_one));
                        stats.max_fn(stats.num_objects+kk,jj) = false_neg(kk,to_consider(which_one));
                        stats.max_inters(stats.num_objects+kk,jj) = inters(kk,to_consider(which_one));
                    end
                end
            end
            stats.obj_classes = [stats.obj_classes; obj_classes(1:size(inters,1))]; 
            stats.num_objects = stats.num_objects + size(jaccards,1);
        else
            error([res_file ' not found']);                
        end
    end
   
    % Check
    if size(stats.obj_classes,1)~=stats.num_objects 
        error('Something went wrong')
    end
   
    stats.mean_n_masks   = mean(stats.all_n_masks);
    
    % ----- Jaccard at instance level (J_i) ----
    % It is the mean best jaccard for all objects
    stats.jaccard_object = mean(stats.max_J);
   
    
    % ----- Compute jaccard at class level (J_c) ----
    class_ids = unique(stats.obj_classes);
    for ii=1:length(class_ids)
        curr_class = class_ids(ii);
        stats.per_class_results{ii}.num_objects = sum(stats.obj_classes==curr_class);
        
        stats.per_class_results{ii}.max_fp     = stats.max_fp(logical(stats.obj_classes==curr_class),:);
        stats.per_class_results{ii}.max_fn     = stats.max_fn(logical(stats.obj_classes==curr_class),:);
        stats.per_class_results{ii}.max_inters = stats.max_inters(logical(stats.obj_classes==curr_class),:);
        stats.per_class_results{ii}.max_J      = stats.max_J(logical(stats.obj_classes==curr_class),:);
        stats.per_class_results{ii}.meanmax    = mean(stats.per_class_results{ii}.max_J);

        stats.per_class_results{ii}.global_fp      = sum(stats.per_class_results{ii}.max_fp,1);
        stats.per_class_results{ii}.global_fn      = sum(stats.per_class_results{ii}.max_fn,1);
        stats.per_class_results{ii}.global_inters  = sum(stats.per_class_results{ii}.max_inters,1);
        
        % Compute per-class total inters, fp, fn
        stats.per_class_results{ii}.global_J = ...
             stats.per_class_results{ii}.global_inters ./...
            (stats.per_class_results{ii}.global_inters+...
             stats.per_class_results{ii}.global_fp    +...
             stats.per_class_results{ii}.global_fn);
    end
    
    % Compute global mean on all classes
    tmp = [];
    for ii=1:length(class_ids)
        tmp = [tmp; stats.per_class_results{ii}.global_J]; %#ok<AGROW>
    end
    stats.jaccard_class = mean(tmp,1);

    save(stats_file,'stats');
    
    % Remove temporal results folder
    rmdir(res_dir,'s');
end
end
