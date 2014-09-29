function a = fitparab(z,ra,rb,theta,filt)
% function [a,b,c] = fitparab(z,ra,rb,theta)
%
% Fit cylindrical parabolas to elliptical patches of z at each
% pixel.  
%
% INPUT
%	z	Values to fit.
%	ra,rb	Radius of elliptical neighborhood, ra=major axis.
%	theta	Orientation of fit (i.e. of minor axis).
%
% OUTPUT
%	a,b,c	Coefficients of fit: a + bx + cx^2
%


% compute the interior quickly with convolutions
a = conv2(z,filt(:,:,1),'same');
%fix border with mex file
a = savgol_border(a, z, ra, rb, theta);

