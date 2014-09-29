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
function [ucm2, ucms, elapsed_time] = img2ucms(I, model, scales, all_parameters, param_multi, align_thr)
% Multiscale hierarchical segmentation
% Pablo Arbelaez 
% arbelaez@berkeley.edu

if nargin<3,
    % This defines the order to project superpixels, so from coarse to fine
    scales = [2, 1, 0.5]; 
end
if nargin<4,
    % BSDS trainval parameters
    all_parameters = [10 60 6 2 0.12]; 
end
if nargin<5,
    param_multi.weights      = ones(1, numel(scales));
    param_multi.thr          = 0.35;
    param_multi.fq           = 8;
end
if nargin<6,
    align_thr = 0.1;
end

param.mult_Pb   = all_parameters(1);
param.sat_sPb   = all_parameters(2);
param.nvec      = all_parameters(3);
param.dthresh   = all_parameters(4);
param.ic_gamma  = all_parameters(5);

% Edit by Jordi Pont-Tuset <jordi.pont@upc.edu>
if size(I,3)==1 % Grayscale image
    I = cat(3,I,I,I); % Otherwise edge_detect fails
end

[param.tx, param.ty, ~] = size(I);

% compute ucms at multiple scales
ucms = cell(numel(scales),1);
elapsed_time = zeros(11,1);
for s = 1:numel(scales),
    param.scale = scales(s);
    [ucms{s}, times] = img2ucm_scale(I, model, param);
    elapsed_time(3*(s-1)+1:3*s)=times;
end

% align ucms
T=tic;
ucms = project_ucms_wrap(ucms, align_thr);
elapsed_time(10)=toc(T);

% combine ucms
T=tic;
ucm2 = ucms2multi(ucms, param_multi);
elapsed_time(11)=toc(T);


%%
function [ucm2, times] = img2ucm_scale(I, model, param)
% compute hierarchical segmentation at each scale independently
I = imresize(I, param.scale, 'lanczos3');
[ucm2, times] = img2ucm(I, model, param.mult_Pb, param.sat_sPb, param.nvec, param.dthresh, param.ic_gamma);

% resample ucm to original image size
T=tic;
if ~isequal([param.tx, param.ty], size(ucm2(3:2:end, 3:2:end))),
    % TODO: Instead of resampling for each threshold, resample leaves and redo mergings
    ucm2 = resample_ucm2(ucm2, [param.tx, param.ty]);
end
times(3)=toc(T);

%%
function [ucm2, times] = img2ucm(I, model, mult_Pb, sat_sPb, nvec, dthresh, ic_gamma)

% edge detection
T=tic;
% TODO: Only once at the beginning 4,2,1 0.5
[E,~,O] = edgesDetect( I, model );  
times(1)=toc(T);

% continuous oriented watershed transform 
T=tic;
[owt2, superpixels] = contours2OWT(E, O);

% globalization
[ sPb_thin] = spectralPb_fast(owt2 * mult_Pb, nvec, ic_gamma, dthresh) / sat_sPb;

% ultrametric contour map with mean pb.
ucm2 = double(ucm_mean_pb( (owt2 + sPb_thin), superpixels) );
times(2)=toc(T);



%%
function ucm2 = ucms2multi(all_ucms, param)

%combine ucms
weights = param.weights;
weights = weights ./ sum(weights);

sz = size(all_ucms);
W_all = repmat(repmat(weights', [1,sz(2)]),[1,1,sz(1)]); W_all = permute(W_all, [3 2 1]);
all_ucms = all_ucms.*W_all;

ucm2_wt = sum(all_ucms,3);

labels2 = bwlabel(ucm2_wt == 0, 8);
labels = labels2(2:2:end, 2:2:end) - 1; % labels begin at 0 in mex file.
ucm2 = double(ucm_mean_pb(ucm2_wt, labels));
bw = (ucm2==0);
ucm2 = apply_sigmoid(ucm2, param.thr, param.fq);
ucm2(bw) = 0;
