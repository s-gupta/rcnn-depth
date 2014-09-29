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
% function masks = cands2labels(candidates,ms)
%
% From a hierarchy (lp: leaves partition and ms: merging sequence) and a
% list of candidates as sets of labels (candidates), we get each candidate
% as a list of labels of the leaves partition lp.
% ------------------------------------------------------------------------
function labels = cands2labels(candidates,ms)
    labels = mex_cands2labels(double(ms),candidates);
end


