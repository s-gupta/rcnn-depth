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

function startup_mcg()
  % Check that 'mcg_root_dir' has been set
  if ~exist(mcg_root_dir,'dir')
      error('Error installing the package, try updating the value of mcg_root_dir in the file "mcg_root_dir.m"')
  end

  % Check that 'mcg_root_dir' has the needed folder
  if ~exist(fullfile(mcg_root_dir,'lib'),'dir')
      error('Error installing the package, the folder "lib" not found, have you compiled it? See mcg_build.m')
  end
  if ~exist(fullfile(mcg_root_dir,'src'),'dir')
      error('Error installing the package, the folder "src" not found')
  end

  % Install own lib
  addpath(mcg_root_dir);
  addpath(fullfile(mcg_root_dir,'lib'));
  addpath(fullfile(mcg_root_dir,'scripts_training'));
  addpath(fullfile(mcg_root_dir,'datasets'));
  addpath(fullfile(mcg_root_dir,'depth_features'));
  addpath(genpath(fullfile(mcg_root_dir,'src')));

  %% Check that the needed functions are compiled
  % Included in our code
  needed_files = {'mex_assess_one_sel','mex_base_perimeters','mex_fast_features',...
                  'mex_fast_intersections', 'mex_fast_reduction',...
                  'mex_get_tree_cands', 'mex_prune_tree_to_regions',...
                  'mex_max_margin', 'mex_hole_filling',...
                  'mex_intersect_hierarchies','mex_cands2masks','mex_cands2labels','mex_ucm2hier',...
                  'mex_eval_masks', 'mex_eval_labels',...
                  'paretofront',...                  % Included from paretofront
                  'mexRF_train', 'mexRF_predict',... % Included from RF_Reg_C
                  'buildW', 'mex_contour_sides', 'ucm_mean_pb',... % Included from BSR
                  %'convConst','gradientMex','imPadMex','imResampleMex','rgbConvertMex',... % Included from piotr_toolbox
                  %'edgesDetectMex',...               % Included from structured_forest
                  };
  for ii=1:length(needed_files)
      if exist(needed_files{ii})~=3 %#ok<EXIST>
          error(['The needed function (' needed_files{ii} ') not found. Have you built the package properly?'])
      end
  end
end
