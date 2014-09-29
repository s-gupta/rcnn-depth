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
function [ucm2_sp] = resample_ucm2_sp(ucm2, superpixels, nthresh)
if nargin<3, nthresh = 100; end

nsp = max(superpixels(:));
mx = max(ucm2(:));
thresh = linspace(mx/nthresh, mx, nthresh)';

ucm2_sp = zeros(size(ucm2));
old_bw = zeros([ (size(ucm2, 1)-1)/2, (size(ucm2, 2)-1)/2]);

for t = 1 : nthresh,
    bw = (ucm2 <= thresh(t) );
    if ~isequal(bw, old_bw),
        labels2 = bwlabel(bw);
        seg = labels2(2:2:end, 2:2:end);
        seg = project_superpixels(seg, superpixels,nsp);
        bdry = seg2bdry(seg);
        old_bw = bw;
    end
    
    ucm2_sp = max(ucm2_sp, thresh(t)*bdry);
end

%%
function seg_sp = project_superpixels(seg, superpixels,nsp)

nsg = max(seg(:));

num1 = nsp + 1;
num2 = nsg + 1;
counts = zeros(num1, num2);

% joint histogram
sumim = 1 + superpixels + seg*num1;
counts(:) = counts(:) + histc(sumim(:), 1:num1*num2);
counts = counts(2:end, 2:end);
[~,id]=max(counts,[],2);
seg_sp = id(superpixels);
