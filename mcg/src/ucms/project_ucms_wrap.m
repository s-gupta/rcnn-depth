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
function ucms_sp = project_ucms_wrap(ucms,thr)

ucms_sp=zeros([size(ucms{1},1),size(ucms{1},2),numel(ucms)]);

for d = 1:numel(ucms),
    ucm_sp = ucms{d};
    for u = d:-1:2
        labels2 = bwlabel(ucms{u-1} <= thr);
        superpixels = labels2(2:2:end, 2:2:end);
        [ucm_sp] = resample_ucm2_sp(ucm_sp, superpixels);
    end
    ucms_sp(:,:,d) = ucm_sp;
end

