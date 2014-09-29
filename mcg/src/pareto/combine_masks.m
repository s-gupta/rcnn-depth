% ------------------------------------------------------------------------ 
%  Copyright (C)
%  Universitat Politecnica de Catalunya BarcelonaTech (UPC) - Spain
% 
%  Jordi Pont-Tuset <jordi.pont@upc.edu>
%  June 2013
% ------------------------------------------------------------------------ 

function [combined_stats,conv_hull_stats,params]= combine_masks(stats1, stats2, measure)
if nargin<3
    measure = 'jaccard_class';
end

n_res1 = length(stats1.mean_n_masks);
n_res2 = length(stats2.mean_n_masks);

%% Combinations
% Single results
for kk=1:n_res1
    combined_stats{kk,n_res2+1} = extract_one_from_all(stats1,kk); %#ok<AGROW>
end
for ii=1:n_res2
    combined_stats{n_res1+1,ii} = extract_one_from_all(stats2,ii); %#ok<AGROW>
end
% Combined results
for kk=1:n_res1
    for ii=1:n_res2
        combined_stats{kk,ii} = join_masks(combined_stats{kk,n_res2+1}, combined_stats{n_res1+1,ii}); %#ok<AGROW>
    end
end
% No results
combined_stats{n_res1+1,n_res2+1}.mean_n_masks = 0;
combined_stats{n_res1+1,n_res2+1}.(measure) = 0;

%% Get convex hull
nn = 1;
all_n_masks = zeros((n_res1+1)*(n_res2+1),1);
xx = zeros((n_res1+1)*(n_res2+1),1);
yy = zeros((n_res1+1)*(n_res2+1),1);
for kk=1:n_res1+1
    for ii=1:n_res2+1
        all_n_masks(nn) = combined_stats{kk,ii}.mean_n_masks;
        xx(nn) = kk;
        yy(nn) = ii;
        nn = nn+1;
    end
end
all_n_masks(end) = [];
xx(end) = [];
yy(end) = [];

nn = 1;
all_J = zeros((n_res1+1)*(n_res2+1),1);
for kk=1:n_res1+1
    for ii=1:n_res2+1
        all_J(nn)   = combined_stats{kk,ii}.(measure);
        nn = nn+1;
    end
end
all_J(end) = [];

if length(all_n_masks)==1
    sel_stats = 1;
else
    % Pareto front
    sel_stats = find(paretoGroup([all_n_masks, -all_J]));
    [~,tmp_ids] = sort(all_n_masks(sel_stats));
    sel_stats = sel_stats(tmp_ids);
    
    % Convex hull to simplify
    sel_stats_hull = convhull(all_n_masks(sel_stats), all_J(sel_stats),'simplify',true);
    sel_stats = sel_stats(sel_stats_hull);
end

for ii=1:length(sel_stats)
    c_hull_stats(ii) = combined_stats{xx(sel_stats(ii)), yy(sel_stats(ii))}; %#ok<AGROW>
    params(ii).params1 = xx(sel_stats(ii)); %#ok<AGROW>
    params(ii).params2 = yy(sel_stats(ii)); %#ok<AGROW>
end


%% Put the convex hull in the same struct
% Sorted in order of n_masks
[~, sel_stats_ids] = sort([c_hull_stats.mean_n_masks]);
params = params(sel_stats_ids);
conv_hull_stats = c_hull_stats(sel_stats_ids(1));
conv_hull_stats.gt_set = stats1.gt_set;
for ii=2:length(sel_stats_ids)
    conv_hull_stats.all_n_masks  = [conv_hull_stats.all_n_masks c_hull_stats(sel_stats_ids(ii)).all_n_masks];
    conv_hull_stats.max_J        = [conv_hull_stats.max_J c_hull_stats(sel_stats_ids(ii)).max_J];
    conv_hull_stats.max_fp       = [conv_hull_stats.max_fp c_hull_stats(sel_stats_ids(ii)).max_fp];
    conv_hull_stats.max_fn       = [conv_hull_stats.max_fn c_hull_stats(sel_stats_ids(ii)).max_fn];
    conv_hull_stats.max_inters   = [conv_hull_stats.max_inters c_hull_stats(sel_stats_ids(ii)).max_inters];
    conv_hull_stats.jaccard_object = [conv_hull_stats.jaccard_object c_hull_stats(sel_stats_ids(ii)).jaccard_object];
    conv_hull_stats.mean_n_masks = [conv_hull_stats.mean_n_masks c_hull_stats(sel_stats_ids(ii)).mean_n_masks];
    conv_hull_stats.jaccard_class  = [conv_hull_stats.jaccard_class c_hull_stats(sel_stats_ids(ii)).jaccard_class];

    for kk=1:length(c_hull_stats(sel_stats_ids(ii)).per_class_results)
        conv_hull_stats.per_class_results{kk}.max_J         = [conv_hull_stats.per_class_results{kk}.max_J c_hull_stats(sel_stats_ids(ii)).per_class_results{kk}.max_J];
        conv_hull_stats.per_class_results{kk}.max_fp        = [conv_hull_stats.per_class_results{kk}.max_fp c_hull_stats(sel_stats_ids(ii)).per_class_results{kk}.max_fp];
        conv_hull_stats.per_class_results{kk}.max_fn        = [conv_hull_stats.per_class_results{kk}.max_fn c_hull_stats(sel_stats_ids(ii)).per_class_results{kk}.max_fn];
        conv_hull_stats.per_class_results{kk}.max_inters    = [conv_hull_stats.per_class_results{kk}.max_inters c_hull_stats(sel_stats_ids(ii)).per_class_results{kk}.max_inters];
        conv_hull_stats.per_class_results{kk}.meanmax       = [conv_hull_stats.per_class_results{kk}.meanmax c_hull_stats(sel_stats_ids(ii)).per_class_results{kk}.meanmax];
  
        conv_hull_stats.per_class_results{kk}.global_fp     = [conv_hull_stats.per_class_results{kk}.global_fp c_hull_stats(sel_stats_ids(ii)).per_class_results{kk}.global_fp];
        conv_hull_stats.per_class_results{kk}.global_fn     = [conv_hull_stats.per_class_results{kk}.global_fn c_hull_stats(sel_stats_ids(ii)).per_class_results{kk}.global_fn];
        conv_hull_stats.per_class_results{kk}.global_inters = [conv_hull_stats.per_class_results{kk}.global_inters c_hull_stats(sel_stats_ids(ii)).per_class_results{kk}.global_inters];
        
        conv_hull_stats.per_class_results{kk}.global_J      = [conv_hull_stats.per_class_results{kk}.global_J c_hull_stats(sel_stats_ids(ii)).per_class_results{kk}.global_J];
    end 
end

%% Some sanity checks
% Compute global mean
tmp = [];
for ii=1:length(conv_hull_stats.per_class_results)
    tmp = [tmp; conv_hull_stats.per_class_results{ii}.global_J]; %#ok<AGROW>
end
assert(isequal(conv_hull_stats.jaccard_class,mean(tmp)));
% Per-class means
for ii=1:length(conv_hull_stats.per_class_results)
    tmp =   conv_hull_stats.per_class_results{ii}.global_inters ./...
            (conv_hull_stats.per_class_results{ii}.global_inters+...
             conv_hull_stats.per_class_results{ii}.global_fp    +...
             conv_hull_stats.per_class_results{ii}.global_fn);
    tmp(isnan(tmp)) = 0;
    assert(isequal(tmp, conv_hull_stats.per_class_results{ii}.global_J))
end
    
%% Debug
% subplot(1,2,1)
% plot(all_n_masks,all_J,'b+')
% hold on 
% plot(all_n_masks(sel_stats),all_J(sel_stats),'k')


end

