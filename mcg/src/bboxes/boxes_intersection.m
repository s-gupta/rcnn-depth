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
function [int_area, int_bbox] = boxes_intersection( bbox1, bbox2 )
    up1    = bbox1(1);
    left1  = bbox1(2);
    down1  = bbox1(3);
    right1 = bbox1(4);

    up2    = bbox2(1);
    left2  = bbox2(2);
    down2  = bbox2(3);
    right2 = bbox2(4);
    
    int_left  = max(left1,left2);
    int_right = min(right1,right2);
    int_up    = max(up1,up2);
    int_down  = min(down1,down2);
    
    if (int_left<=int_right) && (int_up<=int_down)
        int_bbox = [int_up, int_left, int_down, int_right];
        int_area = box_area(int_bbox);
    else
        int_bbox = [];
        int_area = 0;
    end
end

