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
% hierarchy = ucm2hier(ucm)
% ------------------------------------------------------------------------
% Reads a UCM from file or directly as a matrix and converts it to a 
%  hierarcgy of regions (leaves partition and merging_sequence) where
%  we accept non binary mergings, that is, more than one merging for
%  a given threshold.
% ------------------------------------------------------------------------

function hierarchy = ucm2hier(ucm)
% UCM can be a file or the matrix
if ischar(ucm) % ucm refers to a file
    % Get UCM -> Must be saved as a variable named 'ucm2' or 'ucm.strength'
    load(ucm);
    if ~exist('ucm2', 'var')
        ucm2 = ucm.strength;
    end
elseif (size(ucm,1)>2 && size(ucm,2)>2) % It is a full ucm
    ucm2 = ucm;
    clear ucm;
else
    error('UCM type not accepted');
end

% Get leaves segmentation
tmp_ucm = ucm2;
tmp_ucm(1:2:end,1:2:end)=1; % Make the gridbmap connected
labels = bwlabel(tmp_ucm' == 0, 8); % Transposed for the scanning to be from
                                    %   left to right and from up to down
labels = labels';

% ---------------------------

hierarchy.ucm2 = ucm2;
hierarchy.leaves_part = labels(2:2:end, 2:2:end);

% To hierarchy
[hierarchy.ms_matrix, hierarchy.start_ths, hierarchy.end_ths] = mex_ucm2hier(hierarchy.leaves_part, hierarchy.ucm2);

% Store it also as a struct
hierarchy.ms_struct = ms_matrix2struct(hierarchy.ms_matrix);
