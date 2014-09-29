function HHA = saveHHA(imName, C, outDir, D, RD)
% function HHA = saveHHA(imName, C, outDir, D, RD)

% AUTORIGHTS

  if(isempty(D)), D = getImage(imName, 'depth'); end
  if(isempty(RD)), RD = getImage(imName, 'rawdepth'); end
  
  D = double(D)./1000;
  missingMask = RD == 0;
  [pc, N, yDir, h, pcRot, NRot] = processDepthImage(D*100, missingMask, C);
  angl = acosd(min(1,max(-1,sum(bsxfun(@times, N, reshape(yDir, 1, 1, 3)), 3))));
    
  % Making the minimum depth to be 100, to prevent large values for disparity!!!
  pc(:,:,3) = max(pc(:,:,3), 100); 
  I(:,:,1) = 31000./pc(:,:,3); 
  I(:,:,2) = h;
  I(:,:,3) = (angl+128-90); %Keeping some slack
  I = uint8(I);
  
  % Save if can save
  if(~isempty(outDir) && ~isempty(imName)), imwrite(I, fullfile_ext(outDir, imName, 'png')); end
  
  HHA = I;
end
