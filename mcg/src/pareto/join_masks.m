% ------------------------------------------------------------------------ 
%  Copyright (C)
%  Universitat Politecnica de Catalunya BarcelonaTech (UPC) - Spain
% 
%  Jordi Pont-Tuset <jordi.pont@upc.edu>
%  June 2013
% ------------------------------------------------------------------------ 

function [stats, extra] = join_masks(stats1, stats2)
stats.gt_set = stats1.gt_set;

% num_objects
stats.num_objects = stats1.num_objects;
if (stats2.num_objects~=stats1.num_objects)
    error('Incoherent number of number of objects')
end
stats.obj_classes = stats1.obj_classes;

% all_n_masks: [2913x1 double]
if size(stats1.all_n_masks,2)~=size(stats2.all_n_masks,2)
    error('Incoherent number of results')
end
n_results = size(stats2.all_n_masks,2);

stats.all_n_masks = stats1.all_n_masks + stats2.all_n_masks;

stats.max_J = zeros(stats.num_objects, n_results);
stats.max_fp = zeros(stats.num_objects, n_results);
stats.max_fn = zeros(stats.num_objects, n_results);
stats.max_inters = zeros(stats.num_objects, n_results);
extra.which = zeros(stats.num_objects, n_results);
for ii=1:n_results
    [stats.max_J(:,ii), which] = max([stats1.max_J(:,ii), stats2.max_J(:,ii)],[],2);
    extra.which(:,ii) = which;
    for jj=1:stats.num_objects
        if which(jj)==1
            stats.max_fp(jj,ii) = stats1.max_fp(jj,ii);
            stats.max_fn(jj,ii) = stats1.max_fn(jj,ii);
            stats.max_inters(jj,ii) = stats1.max_inters(jj,ii);
        else
            stats.max_fp(jj,ii) = stats2.max_fp(jj,ii);
            stats.max_fn(jj,ii) = stats2.max_fn(jj,ii);
            stats.max_inters(jj,ii) = stats2.max_inters(jj,ii);
        end
    end
end

% Recompute statistics
stats.jaccard_object = mean(stats.max_J);
stats.mean_n_masks   = mean(stats.all_n_masks);


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

    % Compute per-class global J
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