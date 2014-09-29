function [ng1, ng2, dg, sng1, sng2, sdg] = computeDepthCues(z, C, s, depthParam)
% function [ng1, ng2, dg, sng1, sng2, sdg] = computeDepthCues(z, C, s, depthParam)
% Parameters for the local cues
%   C is the camera matrix for the image
%   s is the scaling of the image
% 	depthParam = struct('qzc', 2.9e-5, 'sigmaSpace', 1.40*[1 2 3 4], 'rr', 5*[1 2 3 4], 'sigmaDisparity', [3 3 3 3], 'nori', 8, 'savgolFactor', 1.2);

  % Fill holes if needed
  zf = fillHoles(z, 'bwdist');

  % Convert to point cloud
  [x3 y3 z3] = getPointCloudFromZ(z, C, s);
  [xf3 yf3 zf3] = getPointCloudFromZ(zf, C, s);
  pcf = cat(3, xf3, yf3, zf3);
  pc = cat(3, x3, y3, z3);

  % Compute the local cues
	numScales = length(depthParam.sigmaSpace);
  nori = depthParam.nori;
	
	% Compute the gradients here	
	for i = 1:numScales,
		t = tic; 
		[raw.dg{i} raw.ng{i} raw.sng{i}] = depthCuesHelper(pc, pcf, depthParam.rr(i), depthParam.sigmaSpace(i), depthParam.qzc, depthParam.nori, depthParam.sigmaDisparity(i));
		fprintf('   Time taken for computing ng at scale %d: %0.3f.\n', i, toc(t));
	end

	% Compute the cues and fix their orientations
	for i = 1:numScales,
		ng1{i} = raw.ng{i}; 
		ng1{i}(raw.sng{i} == 1) = 0;
		
		ng2{i} = raw.ng{i}; 
		ng2{i}(raw.sng{i} == -1)= 0;

		%%%%%%%%%   Fix for the channels orientations ....
		ng1{i} = ng1{i}(:,:,mod(nori+1-[1:nori],nori)+1);
		ng2{i} = ng2{i}(:,:,mod(nori+1-[1:nori],nori)+1);
		dg{i} = raw.dg{i}(:,:,mod(nori+1-[1:nori],nori)+1);
		%%%%%%%%%   End for the fix for channel orientations...
	end

  if(nargout() > 3)
    % Smooth the cues
    t = tic();
    gTheta = mod(linspace(pi/2, -pi/2, nori+1), pi);
    gTheta = gTheta(1:end-1);

    sng1 = applySG(ng1, depthParam.rr*depthParam.savgolFactor, gTheta);
    sng2 = applySG(ng2, depthParam.rr*depthParam.savgolFactor, gTheta);
    sdg = applySG(dg, depthParam.rr*depthParam.savgolFactor, gTheta);
    
    sng1 = cat(4,sng1{:});    %NG Convex
    sng2 = cat(4,sng2{:});    %NG Concave
    sdg = cat(4,sdg{:})./100; %DG
    fprintf('   Time taken to apply savgol smoothing %0.2f.\n', toc(t));
  end

  ng1 = cat(4, ng1{:});
  ng2 = cat(4, ng2{:});
  dg = cat(4, dg{:});
end
