function [NR b N err errI nSamples] = computeRegionNormals(pc, superpixels, sp2reg)
%function [NR b N err errI nSamples] = computeRegionNormals(pc, superpixels, sp2reg)
% NR(:,1)'*x + b = 0

  if(~exist('sp2reg','var'))
    sp2reg = eye(max(superpixels(:))) == 1;
  end
  [nSP nReg] = size(sp2reg);
  NR = NaN(3,nReg);
  b = NaN(1,nReg);
  err = NaN(3,nReg);
  
  X = pc(:,:,1);
  for j = 1:nReg,
    regJ = ismember(superpixels, find(sp2reg(:,j)));
    reg = find(regJ & ~isnan(X));
    regJ = find(regJ);
    nSamples(:,j) = [length(reg)./(1+length(regJ)), length(reg)]';
    [NR(:,j) b(j) err(:,j)] = my_fit_plane(pc, reg, 'disparity');
  end

  if(nargin == 2)
    for k = 1:3,
      Ntmp = NR(k,:);
      N(:,:,k) = Ntmp(superpixels);
      errTmp = err(k,:);
      errI(:,:,k) = errTmp(superpixels);
    end
  else
    N = [];
    errI = [];
  end
  %figure(3); imagesc(ismember(superpixels, find(enclosing)));
end
