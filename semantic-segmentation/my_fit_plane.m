function [w b err] = myFitPlane(pc, reg, typ, inlierThresh)
% function [w b err] = myFitPlane(pc, reg, typ, inlierThresh)
% Output: Returns w'x + b = 0;
% Input: pc is Nx3 matrix.

  if(~exist('inlierThresh','var'))
    inlierThresh = 3;
  end
  X = pc(:,:,1);
  Y = pc(:,:,2);
  Z = pc(:,:,3);
  
  A(:,1) = X(reg);
  A(:,2) = Y(reg);
  A(:,3) = Z(reg);
  ind = ~isnan(A(:,1));
  A = A(ind,:);

  if(size(A,1) < 10),
    w = NaN(3,1); b = NaN; err = NaN(3,1);
    return;
  end

  switch typ,
    case 'hack',
        B = A;
        b = 100*ones(size(B,1),1);
        nTmp = linsolve(B,b);
        w = nTmp./norm(nTmp);
        b = -b(1)./norm(nTmp);

    
    case 'disparity',
        % Following XYZ paper which wants to treat disparity as the random variable!
        B = A;
        b = 100*ones(size(B,1),1)./B(:,3);
        B = bsxfun(@times,B,1./B(:,3));
        nTmp = linsolve(B,b);
        w = nTmp./norm(nTmp);
        b = -100./norm(nTmp);

    case 'ransachack',
        B = A(randperm(size(A,1)),:);
        B = A(1:min(500,end),:);
        [wb, P, inliers] = ransacfitplane(B', inlierThresh, false);
        w = wb(1:3)./norm(wb(1:3));
        b = wb(4)./norm(wb(1:3));
  end

  %Reorienting the normals of the floor.
  sg = sign(sum(sign(A*w)));
  if(isnan(sg) || sg == 0)
    sg = 1;
  end
  w = w.*sg;
  b = b*sg;
  
  %Find the errors - fraction of inliers, mean distance from all points, mean distance from inliers.
  dist = A*w+b;
  inliers = find(abs(dist) < inlierThresh);
  err(1) = length(inliers)/size(A,1);
  err(2) = mean(abs(dist));
  err(3) = mean(abs(dist(inliers)));
end
