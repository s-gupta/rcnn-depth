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
function ucm2_sz = resample_ucm2(ucm2, sz)
nthresh=100;

sz2 = 2*sz+1;
thresh = linspace(max(ucm2(:))/nthresh, max(ucm2(:)), nthresh)';

ucm2_sz = zeros(sz2);
old_bw = zeros([ (size(ucm2, 1)-1)/2, (size(ucm2, 2)-1)/2]);

for t = 1 : numel(thresh),
    bw = (ucm2 <= thresh(t) );
    if ~isequal(bw, old_bw),
        labels2 = bwlabel(bw);
        seg = labels2(2:2:end, 2:2:end);
        seg = imresize(seg,sz,'nearest');
        bdry = seg2bdry(seg);
        old_bw = bw;
    end
    ucm2_sz = max(ucm2_sz, thresh(t)*bdry);
end






