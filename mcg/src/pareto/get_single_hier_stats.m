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

function stats = get_single_hier_stats(params,hier_id,gt_set)

if nargin<3
    gt_set = params.gt_set_pareto;
    res_file = params.files.pareto_singles{hier_id};
else
    res_file = strrep(params.files.pareto_singles{hier_id},[params.gt_set_pareto '_pareto_single'],[gt_set '_pareto_single']);
end

% Sampled number of candidates from each hierarchy
n_cands = [10:10:100,200:100:1000,2000:1000:10000];

for ii=1:params.n_r_cand
    n_max_cands(ii) = max(n_cands); %#ok<*AGROW>
end

% Are the results already computed?
if exist(res_file, 'file')
    load(res_file)
    disp(['Loaded: ' res_file '.'])
    recompute = 0;
else
    disp(['RECOMPUTING: ' res_file '.'])
    recompute = 1;
end

if recompute
    % Load which images to consider from the database (train, val, etc.)
    im_ids = database_ids(params.database,gt_set);
    
    for ll=1:params.n_r_cand
        % Store number of regions and gt_set
        stats(ll).n_cands = n_cands;
        stats(ll).gt_set = gt_set;   
        stats(ll).num_objects = 0;   
        stats(ll).obj_classes = [];  
    end
    
    % Sweep all images
    num_images = length(im_ids);
    batches = 20;
    whichOne = imresize(1:batches, [1 num_images], 'nearest');

    for iii = 1:batches,
        indIII = find(whichOne == iii);
        clear Jaccards Inters False_pos False_neg True_areas Obj_classes;

        parfor ii_i = 1:length(indIII),
            ii = indIII(ii_i);
            % Read the UCM as a hierarchy
            hier = ucm2hier(fullfile(params.hier_dirs{hier_id}, [im_ids{ii} '.mat']));
            n_regs = hier.ms_struct(end).parent;
            
            % Get all pairs of neighboring leave regions
            [~, idx_neighbors] = seg2gridbmap(hier.leaves_part);
            K = max(idx_neighbors.matrix_max(:)) + 1;
            neigh_pairs = unique(idx_neighbors.matrix_min+K*idx_neighbors.matrix_max);
            neigh_pairs(neigh_pairs==0) = [];
            neigh_pairs_min = mod(neigh_pairs,K);
            neigh_pairs_max = (neigh_pairs-neigh_pairs_min)/K;

            if isrow(neigh_pairs_min)
                neigh_pairs_min = neigh_pairs_min';
            end
            if isrow(neigh_pairs_max)
                neigh_pairs_max = neigh_pairs_max';
            end

            % Get the 'n_cands' top candidates from each hierarchy
            
            % Singletons
            all_cands = {};
            if n_regs<=n_max_cands(1)
                all_cands{1} = (n_regs:-1:1)';
            else
                all_cands{1} = (n_regs:-1:(n_regs-n_max_cands(1))+1)';
            end

            % Pairs, triplets, etc.
            all_cands(2:params.n_r_cand) = ...
            mex_get_tree_cands(double(hier.leaves_part)-1, double(hier.ms_matrix)-1,...
                                 neigh_pairs_min-1, neigh_pairs_max-1,...
                                 n_max_cands(2:end));
                             
            % Load ground truth
            curr_gt = get_ground_truth( params.database, im_ids{ii} );
            for ll = 1:params.n_r_cand,
                [Jaccards{ii_i}{ll}, Inters{ii_i}{ll}, False_pos{ii_i}{ll}, False_neg{ii_i}{ll}, True_areas{ii_i}{ll}, Obj_classes{ii_i}{ll}] = eval_cands(hier, all_cands{ll}, curr_gt); %#ok<ASGLU>
            end
        end
        fprintf('.');

        for ii_i = 1:length(indIII), 
            ii = indIII(ii_i);
            for ll = 1:params.n_r_cand,
                jaccards = Jaccards{ii_i}{ll};
                inters = Inters{ii_i}{ll};
                false_pos = False_pos{ii_i}{ll};
                false_neg = False_neg{ii_i}{ll};
                true_areas = True_areas{ii_i}{ll};
                obj_classes = Obj_classes{ii_i}{ll};

                % Get best candidate at different number of candidates
                for jj=1:length(stats(ll).n_cands)
                    curr_n_regs = min(stats(ll).n_cands(jj), size(inters,2));
                    to_consider = 1:curr_n_regs;
                    stats(ll).all_n_masks(ii,jj) = length(to_consider); 
                    if (stats(ll).all_n_masks(ii,jj)>0)
                        for kk=1:size(inters,1)
                            [stats(ll).max_J(stats(ll).num_objects+kk,jj), which_one] = max(jaccards(kk,to_consider));
                            stats(ll).max_indicator(stats(ll).num_objects+kk,jj) = to_consider(which_one);
                            stats(ll).max_fp(stats(ll).num_objects+kk,jj)        = false_pos(kk,to_consider(which_one));
                            stats(ll).max_fn(stats(ll).num_objects+kk,jj)        = false_neg(kk,to_consider(which_one));
                            stats(ll).max_inters(stats(ll).num_objects+kk,jj)    = inters(kk,to_consider(which_one));
                        end
                    end
                end
                stats(ll).obj_classes = [stats(ll).obj_classes; obj_classes(1:size(inters,1))'];
                stats(ll).num_objects = stats(ll).num_objects + size(jaccards,1);
            end
        end
    end

    for ll=1:params.n_r_cand,
        stats(ll).jaccard_object = mean(stats(ll).max_J);
        stats(ll).mean_n_masks   = mean(stats(ll).all_n_masks);

        % ----- Compute jaccard at pixel level (J_p) ----
        class_ids = unique(stats(ll).obj_classes);

        % Compute per-class statistics and then mean
        for ii=1:length(class_ids)
            curr_class = class_ids(ii);
            stats(ll).per_class_results{ii}.num_objects = sum(stats(ll).obj_classes==curr_class);

            stats(ll).per_class_results{ii}.max_fp     = stats(ll).max_fp(logical(stats(ll).obj_classes==curr_class),:);
            stats(ll).per_class_results{ii}.max_fn     = stats(ll).max_fn(logical(stats(ll).obj_classes==curr_class),:);
            stats(ll).per_class_results{ii}.max_inters = stats(ll).max_inters(logical(stats(ll).obj_classes==curr_class),:);
            stats(ll).per_class_results{ii}.max_J      = stats(ll).max_J(logical(stats(ll).obj_classes==curr_class),:);
            stats(ll).per_class_results{ii}.meanmax    = mean(stats(ll).per_class_results{ii}.max_J);

            stats(ll).per_class_results{ii}.global_fp      = sum(stats(ll).per_class_results{ii}.max_fp,1);
            stats(ll).per_class_results{ii}.global_fn      = sum(stats(ll).per_class_results{ii}.max_fn,1);
            stats(ll).per_class_results{ii}.global_inters  = sum(stats(ll).per_class_results{ii}.max_inters,1);

            % Compute per-class total inters, fp, fn
            stats(ll).per_class_results{ii}.global_J = ...
                 stats(ll).per_class_results{ii}.global_inters ./...
                (stats(ll).per_class_results{ii}.global_inters+...
                 stats(ll).per_class_results{ii}.global_fp    +...
                 stats(ll).per_class_results{ii}.global_fn);
        end

        % Compute global mean
        tmp = [];
        for ii=1:length(class_ids)
            tmp = [tmp; stats(ll).per_class_results{ii}.global_J];
        end
        stats(ll).jaccard_class = mean(tmp);
    end
    
    % Store
    % fprintf('Done with get_single_hier_stats.\n');
    save(res_file,'stats');
end

