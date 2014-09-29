function fAll = depthFeatures(sp, sp2reg, z, missingMask, C), 
% function fAll = depthFeatures(sp, sp2reg, z, missingMask, C), 
%   sp, sp2reg define the regions over which to compute features
%   z is the depth image in centimeters
%   missingMask boolean variable of whats missing
%   C is the camera matrix

% AUTORIGHTS

  del = 15;
  % z is the depth image as read in ren format, compute the things from it
  % tt = tic();
  [pc, N, yDir, h, pcRot, NRot] = processDepthImage(z, missingMask, C);
  % toc(tt);
  disparity = 1./max(eps, pc(:,:,3));

  % Features come here !!
  % -- mean 
  spArea = accumarray(sp(:), 1);
  regArea = computeSum(sp, sp2reg, ones(size(sp)));
  meanX = computeSum(sp, sp2reg, pcRot(:,:,1))./regArea;
  meanY = computeSum(sp, sp2reg, pcRot(:,:,2))./regArea;
  meanZ = computeSum(sp, sp2reg, pcRot(:,:,3))./regArea;
  meanDisparity = computeSum(sp, sp2reg, disparity)./regArea;

  % -- standard deviation
  stdX = sqrt(computeSum(sp, sp2reg, pcRot(:,:,1).^2)./regArea - meanX.^2);
  stdY = sqrt(computeSum(sp, sp2reg, pcRot(:,:,2).^2)./regArea - meanY.^2);
  stdZ = sqrt(computeSum(sp, sp2reg, pcRot(:,:,3).^2)./regArea - meanZ.^2);
  stdDisparity = sqrt(computeSum(sp, sp2reg, disparity.^2)./regArea - meanDisparity.^2);

  % -- coordinate difference
  xdiff = computeMax(sp, sp2reg, pcRot(:,:,1)) - (-computeMax(sp, sp2reg, -pcRot(:,:,1)));
  ydiff = computeMax(sp, sp2reg, pcRot(:,:,2)) - (-computeMax(sp, sp2reg, -pcRot(:,:,2)));
  zdiff = computeMax(sp, sp2reg, pcRot(:,:,3)) - (-computeMax(sp, sp2reg, -pcRot(:,:,3)));

  % -- hmin, hmax, hdiff, meanH, stdH
  hmin = -computeMax(sp, sp2reg, -h);
  hmax = computeMax(sp, sp2reg, h);
  hdiff = hmax - hmin;
  meanH = computeSum(sp, sp2reg, h)./regArea;
  stdH = sqrt(computeSum(sp, sp2reg, h.^2)./regArea - meanH.^2);

  % -- radial deviation in xz plane
  stdRxz = sqrt(stdX.^2+stdZ.^2);
  stdRxzy = sqrt(stdX.^2+stdZ.^2+stdY.^2);

  angl = acosd(min(1,max(-1,sum(bsxfun(@times, N, reshape(yDir, 1, 1, 3)), 3))));
  
  % -- %age horizontal, %age vertical?
  fracFU = computeSum(sp, sp2reg, double(angl < 0+del))./regArea;
  fracVert = computeSum(sp, sp2reg, double(angl > 90-del & angl < 90+del))./regArea;
  fracFD = computeSum(sp, sp2reg, double(angl > 180-del))./regArea;

  % -- mean angle with gravity
  meanA = computeSum(sp, sp2reg, angl)./regArea;
  stdA = sqrt(computeSum(sp, sp2reg, angl.^2)./regArea - meanA.^2);
  
  % -- max and min deviation in x-z plane, TODO
  % acos^2(t) + bsin^2(t) +2csin(t)cos(t)
  a = stdX.^2; b = stdZ.^2;
  stdXZ = sqrt(computeSum(sp, sp2reg, pcRot(:,:,1).*pcRot(:,:,3))./regArea - meanX.*meanZ);
  c = -stdXZ.^2;
  maxSxz = (a+b)/2 + sqrt(max(0, c.^2 + ((b-a)./2).^2));
  minSxz = (a+b)/2 - sqrt(max(0, c.^2 + ((b-a)./2).^2));
  e = minSxz./max(eps, maxSxz);
  ar = minSxz.*maxSxz;

  fAll = cat(1, regArea, meanX, meanY, meanZ, stdX, stdY, stdZ, xdiff, ydiff, zdiff, hmin, hmax, hdiff, meanH, stdH, stdRxz, stdRxzy, fracFU, fracVert, fracFD, meanA, stdA, maxSxz, minSxz, stdXZ, e, ar, meanDisparity, stdDisparity); 
end

function M = computeSum(sp, sp2reg, m)
% Returns the sum of elements in m for each region
  % tt = tic();
  nSP = max(sp(:));
  ms = accumarray(sp(:), m(:), [nSP 1], @sum, 0);
  M = ms'*sp2reg;
  % fprintf('Time in compute sum %0.3f\n', toc(tt));
end

function M = computeMax(sp, sp2reg, m)
  % tt = tic();
  nSP = max(sp(:));
  ms = accumarray(sp(:), m(:), [nSP 1], @max, -inf);
  ms = bsxfun(@times, ms, sp2reg);
  ms(sp2reg == 0) = -Inf;
  M = max(ms, [], 1);
  % fprintf('Time in compute max %0.3f\n', toc(tt));
end
