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
function sPb_thin = spectralPb_fast(ws_wt2, nvec, ic_gamma, dthresh)
% function sPb_thin = spectralPb_fast(ws_wt2, nvec, ic_gamma, dthresh)
%
% description:
%   fast spectral gradient contours
%
% Jon Barron and Pablo Arbelaez 
% <arbelaez@berkeley.edu>
% Jan 2014

if nargin<4, dthresh = 2; end
if nargin<3, ic_gamma = 0.12; end
if nargin<2, nvec = 6; end

[tx2, ty2] = size(ws_wt2);
tx=(tx2-1)/2; ty=(ty2-1)/2;

l{1} = ws_wt2(1:2:end,2:2:end);
l{2}= ws_wt2(2:2:end,1:2:end);

% build the pairwise affinity matrix
[val,I,J] = buildW(l{1},l{2}, dthresh, ic_gamma);
W = sparse(val,I,J);

[EigVect, EVal] =  ncuts_downsample3(W, nvec, 2, 2, [ty, tx]); 

clear D W opts;

EigVal = diag(EVal);
clear Eval;

EigVal(1:end) = EigVal(end:-1:1);
EigVect(:, 1:end) = EigVect(:, end:-1:1);

vect = zeros(tx, ty, nvec);
for v = 2 : nvec,
    vect(:, :, v) = reshape(EigVect(:, v), [ty tx])';
end
clear EigVect;

%% spectral Pb
for v=2:nvec,
    vect(:,:,v)=(vect(:,:,v)-min(min(vect(:,:,v))))/(max(max(vect(:,:,v)))-min(min(vect(:,:,v))));
end

sPb_thin = zeros(2*tx+1, 2*ty+1);
for v = 1 : nvec
    if EigVal(v) > 0,
        vec = vect(:,:,v)/sqrt(EigVal(v));
        sPb_thin = sPb_thin + seg2bdry_wt(vec, 'doubleSize');
     end
end
sPb_thin = sPb_thin.^(1/sqrt(2));

%%

function [EV, EVal] = ncuts_downsample3(A, NVEC, N_DOWNSAMPLE, DECIMATE, SZ)
% A = affinity matrix
% NEVC = number of eigenvectors (set to 16?)
% N_DOWNSAMPLE = number of downsampling operations (2 seems okay)
% DECIMATE = amount of decimation for each downsampling operation (set to 2)
% SZ = size of the image corresponding to A

A_down = A;
SZ_down = SZ;

Bs = cell(N_DOWNSAMPLE,1);
for di = 1:N_DOWNSAMPLE
    
    % Create a binary array of the pixels that will remain after decimating
    % every other row and column
    [i,j] = ind2sub(SZ_down, 1:size(A_down,1)); 
    do_keep = (mod(i, DECIMATE) == 0) & (mod(j, DECIMATE) == 0);
    
    % Downsample the affinity matrix
    A_sub = A_down(:,do_keep)';
    
    % Normalize the downsampled affinity matrix
    d = (sum(A_sub,1) + eps);
    B = bsxfun(@rdivide, A_sub, d)';
   
    % "Square" the affinity matrix, while downsampling
    A_down = A_sub*B;
    SZ_down = floor(SZ_down / 2);
   
    % Hold onto the normalized affinity matrix for bookkeeping
    Bs{di} = B;
end

% Get the eigenvectors of the Laplacian
%EV = ncuts(A_down, NVEC);
[EV, EVal] = ncuts2(A_down, NVEC);

% diag(EVal)

% Upsample the eigenvectors
for di = N_DOWNSAMPLE:-1:1
    EV = Bs{di} * EV;
end

% whiten the eigenvectors, as they can get scaled weirdly during upsampling
EV = whiten(EV, 1, 0);

function [EV, EVal] = ncuts2(A, n_ev)
[wx, wy] = size(A);
x = 1 : wx;
S = full(sum(A, 1));
D = sparse(x, x, S, wx, wy);
clear S x;


opts.issym=1;
opts.isreal = 1;
opts.disp=0;
[EV, EVal] = eigs((D - A) + (10^-10) * speye(size(D)), D, n_ev, 'sm', opts);
clear D A opts;

v = diag(EVal);
[sorted, sortidx] = sort(v, 'descend');
EV = EV(:,sortidx);
EVal = diag(sorted);

