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
function [ bdry ]  = seg2bdry(seg, fmt)
if nargin<2, fmt = 'doubleSize'; end;

if ~strcmp(fmt,'imageSize') && ~strcmp(fmt,'doubleSize'),
    error('possible values for fmt are: imageSize and doubleSize');
end

[tx, ty, nch] = size(seg);

if nch ~=1, 
    error('seg must be a scalar image');
end

bdry = zeros(2*tx+1, 2*ty+1);

edgels_v = ( seg(1:end-1, :) ~= seg(2:end, :) );
edgels_v(end+1, :) = 0;
edgels_h = ( seg(:, 1:end-1) ~= seg(:, 2:end) );
edgels_h(:, end+1) = 0;

bdry(3:2:end, 2:2:end) = edgels_v;
bdry(2:2:end, 3:2:end) = edgels_h;
bdry(3:2:end-1, 3:2:end-1)= max ( max(edgels_h(1:end-1, 1:end-1), edgels_h(2:end, 1:end-1)), max(edgels_v(1:end-1,1:end-1), edgels_v(1:end-1,2:end)) );

bdry(1, :) = bdry(2, :);
bdry(:, 1) = bdry(:, 2);
bdry(end, :) = bdry(end-1, :);
bdry(:, end) = bdry(:, end-1);

if strcmp(fmt,'imageSize'),
    bdry = bdry(3:2:end, 3:2:end);
end