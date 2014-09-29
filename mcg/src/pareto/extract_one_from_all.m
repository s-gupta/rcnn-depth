% ------------------------------------------------------------------------ 
%  Copyright (C)
%  Universitat Politecnica de Catalunya BarcelonaTech (UPC) - Spain
% 
%  Jordi Pont-Tuset <jordi.pont@upc.edu>
%  June 2013
% ------------------------------------------------------------------------ 

function [ stats ] = extract_one_from_all( stats_all, sel_id)

if sel_id>length(stats_all.mean_n_masks)
    error(['Result ' num2str(sel_id) ' not found']) 
end


stats.gt_set = stats_all.gt_set;
stats.num_objects = stats_all.num_objects;
stats.obj_classes = stats_all.obj_classes;

stats.all_n_masks = stats_all.all_n_masks(:,sel_id);
stats.max_J = stats_all.max_J(:,sel_id);
stats.max_fp = stats_all.max_fp(:,sel_id);
stats.max_fn = stats_all.max_fn(:,sel_id);
stats.max_inters = stats_all.max_inters(:,sel_id);

stats.jaccard_object = stats_all.jaccard_object(:,sel_id);
stats.mean_n_masks = stats_all.mean_n_masks(:,sel_id);


class_ids = unique(stats.obj_classes);

% Compute per-class statistics and then mean
for ii=1:length(class_ids)
    curr_class = class_ids(ii);
    stats.per_class_results{ii}.num_objects = sum(stats.obj_classes==curr_class);

    stats.per_class_results{ii}.max_fp     = stats.max_fp(logical(stats.obj_classes==curr_class),:);
    stats.per_class_results{ii}.max_fn     = stats.max_fn(logical(stats.obj_classes==curr_class),:);
    stats.per_class_results{ii}.max_inters = stats.max_inters(logical(stats.obj_classes==curr_class),:);
    stats.per_class_results{ii}.max_J      = stats.max_J(logical(stats.obj_classes==curr_class),:);
    stats.per_class_results{ii}.meanmax    = mean(stats.per_class_results{ii}.max_J);

    stats.per_class_results{ii}.global_fp      = sum(stats.per_class_results{ii}.max_fp);
    stats.per_class_results{ii}.global_fn      = sum(stats.per_class_results{ii}.max_fn);
    stats.per_class_results{ii}.global_inters  = sum(stats.per_class_results{ii}.max_inters);

    % Compute per-class total inters, fp, fn
    stats.per_class_results{ii}.global_J = ...
         stats.per_class_results{ii}.global_inters ./...
        (stats.per_class_results{ii}.global_inters+...
         stats.per_class_results{ii}.global_fp    +...
         stats.per_class_results{ii}.global_fn);
     
     % NaN means mask has no pixels --> P=1, R=0, F=0
     nan_mark = isnan(stats.per_class_results{ii}.global_J);
     stats.per_class_results{ii}.global_J(nan_mark) = 0;     
end

% Compute global mean
tmp = [];
for ii=1:length(class_ids)
    tmp = [tmp; stats.per_class_results{ii}.global_J]; %#ok<AGROW>
end
stats.jaccard_class = mean(tmp);

end

