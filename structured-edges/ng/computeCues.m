function [dg, ng1, ng2, sdg, sng1, sng2] = computeCues(z, C, s, depthParam)
% Input:
% z is the depth image at cropped size - in units of cm from the camera 
% C is the camera
% s is the scale that you want to use on the image
 
    % Compute the local cues
  [ng1, ng2, dg, raw] = computeDepthCues(pc, pcf, depthParam);

  % Do cue smoothing if needed
  if(nargout() > 3)
    % Do cue smoothing here
  end
end
