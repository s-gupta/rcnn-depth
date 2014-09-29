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
% seg = gridbmap2seg(gridbmap)
%
% From a binary boundary map in the countour grid (as ucm2), get
% the represented partition as matrix of labels in the image plane
%
% INPUT
%	- gridbmap : Binary boundary map on contour grid.
%
% OUTPUT
%	- seg      : Segments labeled from 1..k.
function seg = gridbmap2seg(gridbmap)
tmp = gridbmap;
tmp(1:2:end,1:2:end)=1;
filled_contours = bwlabel(1-tmp',4);
seg = filled_contours(2:2:end,2:2:end)';
