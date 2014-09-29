function [ng1, ng2, dg] = normalCues(z, C, s)
  % Convert to cm to pass into the computeCues function
  z = double(z);
  ss = 1.4;
  depthParam = struct('qzc', 2.9e-5, 'sigmaSpace', ss*[1 2], 'rr', ss*3*[1 2], 'sigmaDisparity', [3 3], 'nori', 4, 'savgolFactor', 1.0);
  % [ng1, ng2, dg, sng1, sng2, sdg] = computeDepthCues(z, C, s, depthParam);
  [ng1, ng2, dg] = computeDepthCues(z, C, s, depthParam);
end
