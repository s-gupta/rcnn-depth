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
%
% This function builds all the MEX files needed.
% Dependencies needed to build: Boost C++ libraries (http://www.boost.org)
%
% ------------------------------------------------------------------------
function mcg_build()
% Check that 'mcg_root_dir' has been set
if ~exist(mcg_root_dir,'dir')
    error('Error building MCG, try updating the value of mcg_root_dir in the file "mcg_root_dir.m"')
end

%% Include the generic paths and files to compile
include{1} = fullfile(mcg_root_dir, 'src', 'aux');  % To get matlab_multiarray.hpp
if (strcmp(computer(),'PCWIN64') || strcmp(computer(),'PCWIN32'))
    include{2} = 'C:\Program Files\boost_1_55_0';  % Boost libraries (change it if necessary)
else
    include{2} = '/opt/local/include/';  % Boost libraries (change it if necessary)
end
include{3} = fullfile(mcg_root_dir, 'src', 'external','piotr_toolbox'); % To build Piotr toolbox

include_str = '';
for ii=1:length(include)
    include_str = [include_str ' -I''' include{ii} '''']; %#ok<AGROW>
end

build_file{1}     = fullfile(mcg_root_dir, 'src', 'cands'    ,'mex_assess_one_sel.cpp');
build_file{end+1} = fullfile(mcg_root_dir, 'src', 'cands'    ,'mex_base_perimeters.cpp');
build_file{end+1} = fullfile(mcg_root_dir, 'src', 'cands'    ,'mex_fast_features.cpp');
build_file{end+1} = fullfile(mcg_root_dir, 'src', 'cands'    ,'mex_fast_intersections.cpp');
build_file{end+1} = fullfile(mcg_root_dir, 'src', 'cands'    ,'mex_fast_reduction.cpp');
build_file{end+1} = fullfile(mcg_root_dir, 'src', 'cands'    ,'mex_get_tree_cands.cpp');
build_file{end+1} = fullfile(mcg_root_dir, 'src', 'cands'    ,'mex_prune_tree_to_regions.cpp');
build_file{end+1} = fullfile(mcg_root_dir, 'src', 'cands'    ,'mex_max_margin.cpp');
build_file{end+1} = fullfile(mcg_root_dir, 'src', 'cands'    ,'mex_hole_filling.cpp');
build_file{end+1} = fullfile(mcg_root_dir, 'src', 'aux'      ,'mex_intersect_hierarchies.cpp');
build_file{end+1} = fullfile(mcg_root_dir, 'src', 'aux'      ,'mex_ucm2hier.cpp');
build_file{end+1} = fullfile(mcg_root_dir, 'src', 'aux'      ,'mex_cands2masks.cpp');
build_file{end+1} = fullfile(mcg_root_dir, 'src', 'aux'      ,'mex_cands2labels.cpp');
build_file{end+1} = fullfile(mcg_root_dir, 'src', 'benchmark','mex_eval_masks.cpp');
build_file{end+1} = fullfile(mcg_root_dir, 'src', 'benchmark','mex_eval_labels.cpp');
build_file{end+1} = fullfile(mcg_root_dir, 'src', 'external' ,'paretofront','paretofront.cpp');

%% Build everything
if ~exist(fullfile(mcg_root_dir, 'lib'),'dir')
    mkdir(fullfile(mcg_root_dir, 'lib'))
end
            
for ii=1:length(build_file)
    eval(['mex ''' build_file{ii} ''' -outdir ''' fullfile(mcg_root_dir, 'lib') '''' include_str])
end

%% Build random forest files
file1   = fullfile(mcg_root_dir, 'src', 'external', 'RF_Reg_C', 'src', 'mex_regressionRF_train.cpp');
file2   = fullfile(mcg_root_dir, 'src', 'external', 'RF_Reg_C', 'src', 'mex_regressionRF_predict.cpp');
dep1    = fullfile(mcg_root_dir, 'src', 'external', 'RF_Reg_C', 'src', 'cokus.cpp');
dep2    = fullfile(mcg_root_dir, 'src', 'external', 'RF_Reg_C', 'src', 'reg_RF.cpp');
o_file1 = fullfile(mcg_root_dir, 'lib', 'mexRF_train');
o_file2 = fullfile(mcg_root_dir, 'lib', 'mexRF_predict');

eval(['mex ' file1 ' ' dep1 ' ' dep2 ' -output ' o_file1 ' -DMATLAB -O'])
eval(['mex ' file2 ' ' dep1 ' ' dep2 ' -output ' o_file2 ' -DMATLAB -O'])

%% Build structured forest files
% eval(['mex ' fullfile(mcg_root_dir, 'src', 'external','structured_forest', 'edgesDetectMex.cpp') ' -outdir ' fullfile(mcg_root_dir, 'lib') include_str])

%% Build piotr_toolbox files
% eval(['mex ' fullfile(mcg_root_dir, 'src', 'external','piotr_toolbox',     'convConst.cpp') ' -outdir ' fullfile(mcg_root_dir, 'lib') include_str])
% eval(['mex ' fullfile(mcg_root_dir, 'src', 'external','piotr_toolbox',   'gradientMex.cpp') ' -outdir ' fullfile(mcg_root_dir, 'lib') include_str])
% eval(['mex ' fullfile(mcg_root_dir, 'src', 'external','piotr_toolbox',      'imPadMex.cpp') ' -outdir ' fullfile(mcg_root_dir, 'lib') include_str])
% eval(['mex ' fullfile(mcg_root_dir, 'src', 'external','piotr_toolbox', 'imResampleMex.cpp') ' -outdir ' fullfile(mcg_root_dir, 'lib') include_str])
% eval(['mex ' fullfile(mcg_root_dir, 'src', 'external','piotr_toolbox', 'rgbConvertMex.cpp') ' -outdir ' fullfile(mcg_root_dir, 'lib') include_str])


%% Build BSR-related files
% 'ucm_mean_pb'
eval(['mex ' fullfile(mcg_root_dir, 'src', 'bsr', 'ucm_mean_pb.cpp') ' -outdir ' fullfile(mcg_root_dir, 'lib')])

% 'buildW'
eval(['mex ' fullfile(mcg_root_dir, 'src', 'bsr', 'buildW.cpp') ' -outdir ' fullfile(mcg_root_dir, 'lib'),...
            ' -I' fullfile(mcg_root_dir,'src','external','BSR','buildW') ' -I' fullfile(mcg_root_dir,'src','external','BSR','buildW','util'),...
            '   ' fullfile(mcg_root_dir,'src','external','BSR','buildW','smatrix.cc'),...
            '   ' fullfile(mcg_root_dir,'src','external','BSR','buildW','ic.cc'),...
            '   ' fullfile(mcg_root_dir,'src','external','BSR','buildW','affinity.cc'),...
            ])
    
% 'mex_contour_sides'
eval(['mex ' fullfile(mcg_root_dir, 'src', 'bsr', 'mex_contour_sides.cpp') ' -outdir ' fullfile(mcg_root_dir, 'lib'),...
            ' -I' fullfile(mcg_root_dir,'src','external','BSR','include'),...
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','concurrent','threads','child_thread.cc'),...
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','concurrent','threads','runnable.cc'),...
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','concurrent','threads','thread.cc'),...
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','concurrent','threads','synchronization','synchronizables','synchronizable.cc'),...
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','concurrent','threads','synchronization','synchronizables','unsynchronized.cc'),...
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','lang','exceptions','ex_bad_cast.cc'),...
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','lang','exceptions','ex_not_found.cc'),...
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','lang','exceptions','ex_not_implemented.cc'),...
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','lang','exceptions','ex_index_out_of_bounds.cc'),...
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','lang','exceptions','ex_invalid_argument.cc'),...
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','lang','exceptions','ex_null_pointer_dereference.cc'),...
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','lang','exceptions','exception.cc'),...
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','lang','exceptions','throwable.cc'),...
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','lang','array.cc'),...                
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','math','random','generators','rand_gen_uniform.cc'),...     
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','math','random','sources','rand_source_default.cc'),...     
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','math','random','sources','rand_source.cc'),...     
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','math','random','sources','mersenne_twister_64.cc'),...     
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','math','random','sources','rand_source_64.cc'),...     
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','math','random','sources','system_entropy.cc'),...     
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','math','random','util','randperm.cc'),...     
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','math','matrices','matrix.cc'),...     
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','math','matrices','exceptions','ex_matrix_dimension_mismatch.cc'),...     
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','math','libraries','lib_image.cc'),...
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','math','libraries','lib_signal.cc'),...                     
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','math','math.cc'),...     
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','math','exact.cc'),...     
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','math','geometry','point_2D.cc'),...     
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','math','geometry','seg_intersect.cc'),...     
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','math','geometry','triangulation.cc'),...  
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','math','geometry','triangle_2D.cc'),...  
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','mlearning','clustering','clusterers','abstract','clusterer.cc'),...  
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','mlearning','clustering','clusterers','abstract','weighted_clusterer.cc'),...  
            '   ' fullfile(mcg_root_dir,'src','external','BSR','src','mlearning','clustering','clusterers','kmeans','basic_clusterer.cc'),...  
            ]);

%% Clear variables
clear build_file file1 file2 dep1 dep2 o_file1 o_file2 ii include include_str

%% Show message
disp('-- Successful compilation of MCG. Enjoy! --')


