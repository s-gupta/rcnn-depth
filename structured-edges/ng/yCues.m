function [y1 y2 y3 angl1 angl2] = yCues(D, C, s)
% function [y1 y2 y3 angl1 angl2] = yCues(D, C, s)
%   Expects D to be in centimetres
%   The camera matrix from somewhere
%   C = getCameraParam('color', 'nyu2');
%   C = fixCamera(C);

% AUTORIGHTS

  % Compute the normals, wrapperComputeNormals expects input in m
  [N, b, pc] = wrapperComputeNormals(double(D)./100, false(size(D)), 3, s, C); 

  % Compute the direction with gravity
	yDirParam.angleThresh = [45 15];
	yDirParam.iter = [5 5];
	yDirParam.y0 = [0 1 0]';
  yDir = getYDir(N, yDirParam);

  % Compute the height above ground
  y0 = [0 1 0]';
  R = getRMatrix2(y0,yDir); pcRot = rotatePC(pc,R'); yMin = min(linIt(pcRot(:,:,2))); 
	y = fillHoles(pcRot(:,:,2)); y = -y; yMin = min(y(:)); 
  if(yMin > -90) yMin = -130; end; y = y-yMin;
  y1 = pc(:,:,2);
  y2 = pcRot(:,:,2);
  y3 = y;
  
  % Compute the angle with gravity vector
  angl1 = sum(bsxfun(@times, N, reshape(yDirParam.y0, 1, 1, 3)), 3);
  angl2 = sum(bsxfun(@times, N, reshape(yDir, 1, 1, 3)), 3);
end
