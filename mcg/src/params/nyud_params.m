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
function params = nyud_params(varargin)
  params = [];
  % Clobber with overrides passed in as arguments
  for i = 1:2:length(varargin)
    key = varargin{i};
    val = varargin{i+1};
    eval(['params.' key ' = val;']);
  end

  params = cv(params, 'gt_set_test', 'test');
  params = cv(params, 'database', 'nyud40Obj');

  % GT sets for training
  params = cv(params, 'gt_set_pareto', 'trainval');
  params = cv(params, 'gt_set_ranking', 'trainval');

  % String identifying the raw hierarchy files
  params = cv(params, 'root_cache_dir', '.');
  params = cv(params, 'raw_hiers', 'hiers');

  % String identifying the current setup (hiers and n_r_cand)
  % to use at different file names
  params = cv(params, 'hiers_id', 'multi_3sc_u_3r');  % multi and 3 scales 'up'

  % String identifying the current Pareto chosen point
  params = cv(params, 'pareto_id', '20k');

  % Maximum number of regions per candidate
  params = cv(params, 'n_r_cand', 3);

  % To identify whether we want to use depth based features or not
  params = cv(params, 'feature_id', 'no_depth');
  params = cv(params, 'depth_features', false);
  params = cv(params, 'camera_matrix', []);

  % Hierarchies to combine
  params.hiers    = {'multi','scale_2.0', 'scale_1.0', 'scale_0.5'};


  % Level of overlap to erase duplicates
  params = cv(params, 'J_th', 0.95);

  % Number of samples per image that the ranking training will get
  params = cv(params, 'n_samples', 30);

  % Max margin parameter
  params = cv(params, 'theta', 0.7);


  %% Do not edit the lines below
  % Where to find the hierarchies
  for ii=1:length(params.hiers)
      params.hier_dirs{ii} = fullfile(params.root_cache_dir, 'output', 'ucm', params.hiers{ii});
  end
  
  % Cache directory
  params.cacheDir = fullfile(params.root_cache_dir, 'regions', 'cache');
  exists_or_mkdir(params.cacheDir);
  
  % String identifying the current final global setup at MCG
  params.mcg_id = [params.raw_hiers '_' params.hiers_id '_' params.pareto_id '_' params.feature_id];

  % ----- Results files -----
  % Global Pareto file
  params.files.pareto_full = fullfile(params.root_cache_dir, 'regions', 'pareto', [params.raw_hiers '_' params.hiers_id '_nr_' num2str(params.n_r_cand) '_' params.gt_set_pareto '_pareto_full.mat']);

  % Single-scale results
  for ii=1:length(params.hiers)
      params.files.pareto_singles{ii} = fullfile(params.root_cache_dir, 'regions', 'pareto', [params.raw_hiers '_' params.hiers{ii} '_nr_' num2str(params.n_r_cand) '_' params.gt_set_pareto '_pareto_single.mat']);
  end

  % Parameters of chosen Pareto point
  params.files.output_dir = fullfile(params.root_cache_dir, 'output', 'regions', params.mcg_id);
  exists_or_mkdir(params.files.output_dir);
  
  params.files.pareto_point = fullfile(params.root_cache_dir, 'regions', 'pareto', [params.mcg_id '_' params.gt_set_pareto '_pareto_point.mat']);
  exists_or_mkdir(fullfile(params.root_cache_dir, 'regions', 'pareto'));

  % File where full features will be stored
  params.files.features_file = fullfile(params.root_cache_dir, 'regions', 'features', [params.mcg_id '_features_' params.gt_set_ranking  '.mat']);
  exists_or_mkdir(fullfile(params.root_cache_dir, 'regions', 'features'));

  % File where the trained classifier will be stored
  params.files.trained_classifier = fullfile(params.root_cache_dir, 'regions', 'classifiers',[params.mcg_id '_rand_forest_' params.gt_set_ranking  '.mat']);
  exists_or_mkdir(fullfile(params.root_cache_dir, 'regions', 'classifiers'));
end

% -------------------------------------------------------------------
% Does nothing if conf.key exists, otherwise sets conf.key to val
function conf = cv(conf, key, val)
  try
    eval(['conf.' key ';']);
  catch
    eval(['conf.' key ' = val;']);
  end
end
