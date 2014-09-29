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
function rank_training(params, cache_dir)

n_hiers = length(params.hiers);

%% Load Pareto parameters
n_cands = loadvar(params.files.pareto_point,'n_cands');

% Some checks
assert(n_hiers==size(n_cands,2));
assert(params.n_r_cand==size(n_cands,1));

%% Gather all features and quality of the training set candidates
if exist(params.files.features_file,'file')
    load(params.files.features_file);
    disp(['Loaded: ' params.files.features_file '.'])
else
    disp(['RECOMPUTING: ' params.files.features_file '.'])

    % features = [];
    % jaccards = [];

    % Load which images to consider from the params.database (train, val, etc.)
    im_ids = database_ids(params.database,params.gt_set_ranking);

    num_images = length(im_ids);
    parfor im_id = 1:num_images
        out = cache_mcg_features(params, [], im_ids{im_id}, cache_dir);
        f_lp = out.f_lp; f_ms = full(out.f_ms); feats = out.feats;
        bboxes = out.bboxes; red_cands = out.red_cands; 
        b_feats_intersections = full(out.b_feats_intersection); 

        if(params.depth_features)
          tt = tic();
          sp = f_lp;
          sp2regC = cands2labels(red_cands, f_ms);
          sp2reg = false(max(sp(:)), length(sp2regC));
          for i = 1:length(sp2regC), sp2reg(sp2regC{i},i) = true; end
          D = getImage(im_ids{im_id}, 'depth');
          D = double(D)./1000;
          C = params.camera_matrix;
          missingMask = getImage(im_ids{im_id}, 'rawdepth') == 0;
          fdepth = depthFeatures(sp, sp2reg, D, missingMask, C);
          feats = cat(2, feats, fdepth');
          fprintf('%s: Time for depth features: %0.3f\n', im_ids{im_id}, toc(tt));
        end

        % Eval candidates
        gt = get_ground_truth(params.database, im_ids{im_id});
        hier = struct('leaves_part', f_lp, 'ms_struct', ms_matrix2struct(f_ms));
        jacc = eval_cands(hier, red_cands, gt);
        max_jacc = max(jacc,[],1)';
        
        % fprintf('.');
        % Sample candidates
        % Get the optimum candidates for each object
        % to ensure they are on the training set
        [~,ids] = max(jacc,[],2);
        if (length(ids)>params.n_samples)  % More objects than samples asked
            ids = ids(1:params.n_samples);
        end
        ids_rest = setdiff(1:size(jacc,2),ids);
        if (size(jacc,2)-length(ids))>params.n_samples
            ids_rest = ids_rest(randperm(length(ids_rest),params.n_samples-length(ids)));
        end
        if isrow(ids_rest)
            ids_rest = ids_rest';
        end
        sel_ids = [ids; ids_rest];

        % Store
        % features = [features; feats(sel_ids,:)]; %#ok<AGROW>
        % jaccards = [jaccards; max_jacc(sel_ids,:)]; %#ok<AGROW>
        
        features{im_id} = feats(sel_ids,:); 
        jaccards{im_id} = max_jacc(sel_ids,:);
    end
    features = cat(1, features{:});
    jaccards = cat(1, jaccards{:});

    save(params.files.features_file,'features','jaccards');
end

%% Train the random forest
if exist(params.files.trained_classifier,'file')
    disp(['Already trained: ' params.files.trained_classifier '.'])
else
    disp(['TRAINING: ' params.files.trained_classifier '.'])

    % Train and save the result
    rf = regRF_train(features,jaccards,50); %#ok<NASGU>
    disp('Training done')

    % Save the result
    save(params.files.trained_classifier,'rf');
end
