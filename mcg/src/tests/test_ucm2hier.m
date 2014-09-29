% ------------------------------------------------------------------------ 
%  Copyright (C)
%  Universitat Politecnica de Catalunya BarcelonaTech (UPC) - Spain
% 
%  Jordi Pont-Tuset <jordi.pont@upc.edu>
%  June 2013
% ------------------------------------------------------------------------ 

%% Dummy

ucm2 = [ 0   0   0   0   0   0   0
         0   0   0   0   0.3 0   0
         0   0.4 0   0.4 0   0.5 0
         0   0   0   0   0.5 0   0
         0   0   0   0   0   0   0];
     
hier = ucm2hier(ucm2);
assert(isequal(hier.leaves_part,[1 1 2; 3 3 4]))
assert(isequal(hier.ms_matrix,[1 2 5; 3 5 6; 4 6 7]))
assert(isequal(hier.start_ths,[0 0 0 0 0.3 0.4 0.5]'))
assert(isequal(hier.end_ths,[0.3 0.3 0.4 0.5 0.4 0.5 1]'))


ucm2 = [ 0   0   0   0   0   0   0
         0   0   0   0   0.5 0   0
         0   0.4 0   0.4 0   0.5 0
         0   0   0.3 0   0.3 0   0
         0   0   0   0   0   0   0];
     
hier = ucm2hier(ucm2);
assert(isequal(hier.leaves_part,[1 1 2; 3 4 5]))
assert(isequal(hier.ms_matrix,[3 4 5 6; 1 6 0 7; 2 7 0 8]))
assert(isequal(hier.start_ths,[0 0 0 0 0 0.3 0.4 0.5]'))
assert(isequal(hier.end_ths,[0.4 0.5 0.3 0.3 0.3 0.4 0.5 1]'))


ucm2 = [ 0   0   0   0   0   0   0
         0   0   0   0   0   0   0
         0   0   0   0   0   0   0
         0   0   0   0   0   0   0
         0   0.9 0   0.9 0   0.9 0
         0   0   0   0   0.5 0   0
         0   0   0   0   0   0   0];

hier = ucm2hier(ucm2);

assert(isequal(hier.leaves_part,[1 1 1; 1 1 1; 2 2 3]))
assert(isequal(hier.ms_matrix,[2 3 4; 1 4 5]))


ucm2 = [ 0   0   0   0   0   0   0    0  0
         0   0   0.1 0   0.2 0   0.1  0  0
         0   0   0   0   0   0   0    0  0];

hier = ucm2hier(ucm2);

assert(isequal(hier.leaves_part,[1 2 3 4]))
assert(isequal(hier.ms_matrix,[1 2 5; 3 4 6; 5 6 7]))
assert(isequal(hier.start_ths,[0 0 0 0 0.1 0.1 0.2]'))
assert(isequal(hier.end_ths,[0.1 0.1 0.1 0.1 0.2 0.2 1]'))


ucm2 = [ 0   0   0   0   0    0   0    0  0
         0   0   0.1 0   0.05 0   0.1  0  0
         0   0   0   0   0    0   0    0  0];

hier = ucm2hier(ucm2);

assert(isequal(hier.leaves_part,[1 2 3 4]))
assert(isequal(hier.ms_matrix,[2 3 0 5; 1 4 5 6]))
assert(isequal(hier.start_ths,[0 0 0 0 0.05 0.1]'))
assert(isequal(hier.end_ths,[0.1 0.05 0.05 0.1 0.1 1]'))


%% Real
load(fullfile(root_dir,'datasets','pascal2012','gPb_mUCM','multi','2007_000033.mat'));
hier  = ucm2hier(ucm2);
assert(length(unique(hier.leaves_part)) + size(hier.ms_matrix,1)==7543)
assert(isequal(unique(ucm2), unique(hier.start_ths)))

