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
% Create random divisions of the training set for cross-validation
%
% ------------------------------------------------------------------------

% Create random divisions of the training set
id = '3';
percentage_train = 0.5;
full_set = 'train2012';

% Get full ids 
im_ids = database_ids('pascal2012',full_set);

ids_a = sort(randperm(length(im_ids),floor(percentage_train*length(im_ids))));
ids_b = setdiff(1:length(im_ids),ids_a);
im_ids_a = im_ids(ids_a);
im_ids_b = im_ids(ids_b);

% Write to file
file_a = fullfile(root_dir,'datasets', 'pascal2012','gt_sets',[full_set '_' id 'a.txt']);
file_b = fullfile(root_dir,'datasets', 'pascal2012','gt_sets',[full_set '_' id 'b.txt']);

dlmwrite(file_a,im_ids_a,'')
dlmwrite(file_b,im_ids_b,'')