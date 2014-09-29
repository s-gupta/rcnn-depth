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
function params = sf_mUCM_multi_3sc_u_4r_12k_params(  )

% Pascal database
params.database = 'pascal2012';

% GT sets for training
params.gt_set_pareto  = 'train2012';
params.gt_set_ranking = 'train2012';

% GT set for testing
params.gt_set_test = 'val2012';

% String identifying the raw hierarchy files
params.raw_hiers = 'sf_mUCM';

% String identifying the current setup (hiers and n_r_cand)
% to use at different file names
params.hiers_id = 'multi_3sc_u_4r';  % multi and 3 scales 'up'

% String identifying the current Pareto chosen point
params.pareto_id = '12k';

% Maximum number of regions per candidate
params.n_r_cand = 4;

% Hierarchies to combine
params.hiers    = {'multi','scale_0.50','scale_1.00','scale_2.00'};

% Level of overlap to erase duplicates
params.J_th = 0.95;

% Number of samples per image that the ranking training will get
params.n_samples = 30;

% Max margin parameter
params.theta = 0.7;


%% Do not edit the lines below
% Where to find the hierarchies
for ii=1:length(params.hiers)
    params.hier_dirs{ii} = fullfile(root_dir,'datasets',params.database,params.raw_hiers,params.hiers{ii});
end

% String identifying the current final global setup at MCG
params.mcg_id = [params.raw_hiers '_' params.hiers_id '_' params.pareto_id];

% ----- Results files -----
% Global Pareto file
params.files.pareto_full = fullfile(root_dir,'datasets',params.database,'pareto',[params.raw_hiers '_' params.hiers_id '_nr_' num2str(params.n_r_cand) '_' params.gt_set_pareto '_pareto_full.mat']);

% Single-scale results
for ii=1:length(params.hiers)
    params.files.pareto_singles{ii} = fullfile(root_dir,'results',params.database,'pareto',[params.raw_hiers '_' params.hiers{ii} '_nr_' num2str(params.n_r_cand) '_' params.gt_set_pareto '_pareto_single.mat']);
end

% Parameters of chosen Pareto point
params.files.pareto_point = fullfile(root_dir,'datasets',params.database,'pareto',[params.mcg_id '_' params.gt_set_pareto '_pareto_point.mat']);

% File where full features will be stored
params.files.features_file = fullfile(root_dir,'datasets',params.database,'features',[params.mcg_id '_features_' params.gt_set_ranking  '.mat']);

% File where the trained classifier will be stored
params.files.trained_classifier = fullfile(root_dir,'datasets',params.database,'classifiers',[params.mcg_id '_rand_forest_' params.gt_set_ranking  '.mat']);


% ---- Ensure needed folders ----
if ~exist(fullfile(root_dir,'results',params.database,'pareto'),'dir')
    mkdir(fullfile(root_dir,'results',params.database,'pareto'))
end
if ~exist(fullfile(root_dir,'results',params.database,'cands'),'dir')
    mkdir(fullfile(root_dir,'results',params.database,'cands'))
end
if ~exist(fullfile(root_dir,'datasets',params.database,'classifiers'),'dir')
    mkdir(fullfile(root_dir,'datasets',params.database,'classifiers'));
end
if ~exist(fullfile(root_dir,'datasets',params.database,'features'),'dir')
    mkdir(fullfile(root_dir,'datasets',params.database,'features'));
end