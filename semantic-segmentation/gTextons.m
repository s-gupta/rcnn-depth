function [map] = gTextons(pc, N, yDir, param)
% function [] = gTextons(pc, N, yDir, param)
%   param.nbins for each angl and height
% param,yMinMean    [= -130] the value to use for the height above floor if unable to find floor.
% param.yMinOutlier [= -90] the value to determine if the floor found is correct or not.

  y0 = [0 1 0]';
  R = getRMatrix(y0,yDir);
  NRot = rotatePC(N,R');
  pcRot = rotatePC(pc,R');

  angl = acosd(sum(bsxfun(@times, N, reshape(yDir,[1 1 3])),3));
  y = -(pcRot(:,:,2));
  yMin = prctile(y(:),2);
  if(yMin > param.yMinOutlier)
    yMin = param.yMinMean;
  end
  y = y-yMin;
  
  map = toGTextons(y, angl, param.nbins);
end

function st = toGTextons(y, angl, nbins)
% function shapeTextons = toGTextons(y, angl, nbins)
% Input: y is the height above ground in metres
%    angl is the angl with the gravity vector from 0 to 180 degree
%    nbins is the number of bins
% Output: shapeTextons is the output quantized shape textons.
  yQ = ceil(max(min(y,200),1)/200*(nbins-1));
  yQ(isnan(yQ)) = nbins;
  
  anglQ = ceil(max(min(angl,180),1)/180*(nbins-1));
  anglQ(isnan(anglQ)) = nbins;
  
  st = yQ+nbins*(anglQ-1);
end
