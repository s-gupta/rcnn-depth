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
function mask = box2mask( bbox, sx, sy )
%     up    = bbox(1);
%     left  = bbox(2);
%     down  = bbox(3);
%     right = bbox(4);
    
    horiz = zeros(1,sy);
    vert  = zeros(1,sx);
    horiz(bbox(2):bbox(4)) = 1;
    vert(bbox(1):bbox(3)) = 1;
    mask = vert'*horiz;
end

