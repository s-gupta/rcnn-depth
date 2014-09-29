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
function bbox = mask2box( mask )
    tmp = regionprops(double(mask),'BoundingBox'); % Double to force a single bbox
    bbox(1) = tmp.BoundingBox(2)+0.5;
    bbox(2) = tmp.BoundingBox(1)+0.5;
    bbox(3) = tmp.BoundingBox(2)+tmp.BoundingBox(4)-0.5;
    bbox(4) = tmp.BoundingBox(1)+tmp.BoundingBox(3)-0.5;
end

