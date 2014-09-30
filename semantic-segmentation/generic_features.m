function features = genericFeatures(superpixels, sp2reg, param, data)

  %param.selN = 4;
  %param.localNormalExt = 'oPbV5_6';

  [nSP nReg] = size(sp2reg);
  spArea = histc(superpixels(:), 1:nSP)';
  regArea = spArea*sp2reg;

  % load the point cloud.
  pc = data.pc;
  N = data.normals;
  yDir = data.yDir;
  bgMax = data.bgMax;
  bgOriented = data.bgOriented;   %= data.colorCues(:,:,:,1:3);
  zgMax = data.zgMax;       %max(data.depthCues(:,:,:,1:15), [], 3);


  % Rotate the point cloud and the normals
  y0 = [0 1 0]';
  R = getRMatrix(y0,yDir);
  NRot = rotatePC(N,R');
  pcRot = rotatePC(pc,R');

  % Find normals for each region
  [NR b gr grr nSamples] = computeRegionNormals(pc, superpixels, sp2reg);
  NRRot = R'*NR;
  
  %Feature 1 angle with the gravity vector.
  tic;
  angl = acosd(sum(bsxfun(@times, N, reshape(yDir,[1 1 3])),3));
  af = @(i) angleFeatures(angl, NR, yDir, i, sp2reg, superpixels);
  angleF = arrayfun(af, 1:nReg,'UniformOutput', false);
  angleF = cell2mat(angleF);
  F{1} = angleF;
  fprintf('Angle Features : %0.3f seconds.\n',toc);

  %Planarity - variance of the points on the plane.
  tic;
  pf = @(w,b,i) planarityFeatures(pc, w, b, i, sp2reg, superpixels);
  allDist = arrayfun(pf, mat2cell(NR,3,ones(1,size(NR,2))), b, 1:nReg,'UniformOutput', false);
  allDist = cell2mat(allDist);
  F{2} = allDist;
  fprintf('Planarity Features : %0.3f seconds.\n',toc);
  
  %Feature 3 is the normal and the bias for the plane, and the area in the image.
  tic;
  F{3} = [NRRot; b./100; log(regArea)];
  fprintf('Equation of plane features : %0.3f seconds.\n',toc);


  %Feature 4 is the range of X span, Z span and Ymin and Ymax in terms of aligned coordinates.
  tic;
  xR = pcRot(:,:,1);
  yR = -pcRot(:,:,2);
  zR = pcRot(:,:,3);
  
  ef = @(i) extentFeatures(xR, yR, zR, i, sp2reg, superpixels);
  eF = arrayfun(ef, 1:nReg, 'UniformOutput', false);
  F{4} = cell2mat(eF);
  fprintf('Extent Features : %0.3f seconds.\n',toc);


  %area features
  tic;
  x = pc(:,:,1);
  y = -pc(:,:,2);
  z = pc(:,:,3);
  F{5} = areaFeatures(x, y, z, N, yR, yDir, sp2reg, superpixels); 
  fprintf('Area Features : %0.3f seconds.\n',toc);
  %F{6} = areaFeaturesBatch(x, y, z, N, yR, yDir, sp2reg, superpixels); 


  %Contour features.
  tic;
  %load something, right?
  % dt.features = getFeaturesLearnWtV1OnlyMax(im_i);
  % dt.features = reshape(dt.features',[size(z) 28]);
  % zSignal(:,:,1) = mean(dt.features(:,:,1:5),3);
  % zSignal(:,:,2) = mean(dt.features(:,:,6:10),3);
  % zSignal(:,:,3) = mean(dt.features(:,:,11:15),3);
  % zSignal(:,:,4) = mean(dt.features(:,:,16:18),3);
  zSignal(:,:,1) = mean(zgMax(:,:,1:5),3);
  zSignal(:,:,2) = mean(zgMax(:,:,6:10),3);
  zSignal(:,:,3) = mean(zgMax(:,:,11:15),3);
  zSignal(:,:,4) = mean(bgMax,3);
  zGood = ~imerode(isnan(z),strel('disk',10));
  cf = @(i) contourFeatures(zSignal, zGood, i, sp2reg, superpixels);
  contourF = arrayfun(cf, 1:nReg, 'UniformOutput', false);
  F{6} = cell2mat(contourF);
  fprintf('Contour Energy Features : %0.3f seconds.\n',toc);

  %Clipping features
  tic;
  cf = @(i) clippingFeatures(i, sp2reg, superpixels);
  clipF = arrayfun(cf, 1:nReg, 'UniformOutput', false);
  F{7} = cell2mat(clipF);
  fprintf('Is Clipped Features : %0.3f seconds.\n',toc);

  %BG Oriented energy feature
  tic;
  % bg = load(sprintf('%s/alignedImagesF/bgV11/img_%04d_bgV11.mat',c.resultsDir,im_i),'BG');
  oBG = orientedEnergy(bgOriented, sp2reg, superpixels);
  F{8} = oBG;
  fprintf('Oriented Contour Energy Features : %0.3f seconds.\n',toc);

  %Average orientation of patches around me
  tic;
  oC = @(i) orientationContext(angl, i, sp2reg, superpixels);
  oContextF = arrayfun(oC, 1:nReg, 'UniformOutput', false);
  F{9} = cell2mat(oContextF);
  fprintf('Orientation Context Feature : %0.3f seconds.\n',toc);

  features = cat(1,F{:});
end

function f = orientationContext(angl, i, sp2reg, superpixels)
  ind = ismember(superpixels, find(sp2reg(:,i)));
  bb = regionprops(ind, 'boundingbox');
  bb = round(cat(1,bb.BoundingBox));
  bb = [min(bb(:,1)) min(bb(:,2)) max(bb(:,1)+bb(:,3))-min(bb(:,1)), max(bb(:,2)+bb(:,4))-min(bb(:,2))];
  bb = [bb(2) bb(1) bb(4) bb(3)];
  m = zeros(size(superpixels));
  for i = [-1 0 1],
    for j = [-1 0 1],
      m(max(bb(1)+i*bb(3),1):min(bb(1)+(i+1)*bb(3)-1,size(m,1)), max(bb(2)+j*bb(4),1):min(bb(2)+(j+1)*bb(4)-1,size(m,2))) = (j+2)+3*(i+1);
    end
  end
  m = m.*~isnan(angl);
  f = accumarray(m(m>0),angl(m>0),[9 1],@mean,0);
  g = accumarray(m(m>0),angl(m>0),[9 1],@median,0);
  f = [f;g];
end


function f = orientedEnergy(bg, sp2reg, superpixels)
  ep = 0.1;
  spArea = accumarray(superpixels(:),1)';
  regArea = spArea*sp2reg;
  for i = 1:size(bg,3),
    bgi = bg(:,:,i);
    spEnergy = accumarray(superpixels(:),bgi(:));
    regEnergy(i,:) = spEnergy'*sp2reg;
  end
  regEnergy = bsxfun(@times, regEnergy, 1./regArea);
  f = regEnergy;
  g = bsxfun(@times, regEnergy, 1./sqrt(sum(regEnergy.^2,1)+ep^2)); 
  f = [f; g];
end

function f = clippingFeatures(i, sp2reg, superpixels)
  ind = ismember(superpixels, find(sp2reg(:,i))); 
  t = regionprops(ind, 'Perimeter');
  p = sum(cat(1,t.Perimeter))+1;
  tol = 2;
  onEdge = sum(any(ind(1:tol,:),1)) + sum(any(ind(end-tol+1:end,:),1)) + sum(any(ind(:,1:tol),2)) + sum(any(ind(:,end-tol+1:end),2)) + 1;
  f = [p onEdge log(p) log(onEdge) onEdge./p]';
end

function f = contourFeatures(zg, zGood, i, sp2reg, superpixels)
  %zg is a 3 channel thing.
  rad = 5;
  regO = ismember(superpixels,find(sp2reg(:,i)));
  
  %imerode the region
  reg{1} = imerode(regO,strel('disk',rad));
  reg{2} = imdilate(regO,strel('disk',rad));
  %reg{3} = imdilate(regO,strel('disk',rad*2));
  reg{3} = true(size(superpixels));
  for i = 1:length(reg),
    zgi = squeeze(sum(sum(bsxfun(@times, zg, reg{i} & zGood),1),2));
    goodArea = nnz(zGood & reg{i});
    totalArea = nnz(reg{i})+1;
    f{i} = [zgi./totalArea];
  end
  f{4} = f{2}-f{1};
  f{5} = f{3}-f{2};
  %f{7} = f{4}-f{2};
  f = cell2mat(f');
end

function f = angleFeatures(angl, NR, yDir, i, sp2reg, superpixels)
  del = 15;
  %Angle of the whole plane with vertical
  f(1) = acosd(NR(:,i)'*yDir);

  reg = ismember(superpixels, find(sp2reg(:,i))) & ~isnan(angl);
  %Mean Local Angle
  f(2) = mean(angl(reg));

  %Median Local Angle
  f(3) = median(angl(reg));

  % %age horizontal, %age vertical?
  f(4) = nnz(angl(reg) < 0+del)./nnz(reg);
  f(5) = nnz(angl(reg) > 180-del)./nnz(reg);
  f(6) = nnz(angl(reg) > 90-del & angl(reg) < 90+del)./nnz(reg);
  f = f';
end

function f = areaFeatures(x, y, z, N, yR, yDir, sp2reg, superpixels)
  %Calculate the angle of each patch with the viewing direction.
  fx = 520;
  fy = 520;
  del = 15;
  nSP = max(superpixels(:));
  angl = abs(sum(N.*cat(3,x,y,z),3)./sqrt(sum(cat(3,x,y,z).^2,3)));
  %Area of the surface
  dArea = (z.*z./angl)./(fx.*fy);
  ind = find(~isnan(dArea));

  spArea = accumarray(superpixels(:), superpixels(:) > 0,[nSP 1]);
  goodSPArea = accumarray(superpixels(ind), superpixels(ind) > 0,[nSP 1]);
  ar = accumarray(superpixels(ind),dArea(ind), [nSP 1]);
  f(1,:) = ar'*sp2reg; %Surface area of the region in space using only what is available for measrement
  f(2,:) = ar'*sp2reg./(goodSPArea'*sp2reg).*(spArea'*sp2reg); %Surface area of region scaled by the amount of missing pixels in the superpixels.
  
  %Horizontal area
  %Find the horizontal pixels and find there area
  angl = max(min(sum(bsxfun(@times, reshape(yDir,[1 1 3]), N),3),1),-1);
  angl = acosd(angl);

  %Horizontal Area and height of horizontal features.
  ind = angl < 0+del & ~isnan(dArea);
  ar = accumarray(superpixels(ind),dArea(ind), [nSP 1], @sum, 0);
  f(3,:) = ar'*sp2reg; 
  ind = angl < 0+del & ~isnan(yR);
  ar = accumarray(superpixels(ind),yR(ind),[nSP 1],@sum, 0);
  cnt = accumarray(superpixels(ind),superpixels(ind) > 0,[nSP 1],@sum, 0);
  f(4,:) = (ar'*sp2reg)./(max(cnt'*sp2reg,1));
  for i = 1:size(sp2reg,2),
    reg = ismember(superpixels, find(sp2reg(:,i))) & ind;
    f(5,i) = median([yR(reg); 0]);
  end
  
  ind = angl > 180-del & ~isnan(dArea);
  ar = accumarray(superpixels(ind),dArea(ind), [nSP 1], @sum, 0);
  f(6,:) = ar'*sp2reg; 
  ind = angl > 180-del & ~isnan(yR);
  ar = accumarray(superpixels(ind),yR(ind),[nSP 1],@sum, 0);
  cnt = accumarray(superpixels(ind),superpixels(ind) > 0,[nSP 1],@sum, 0);
  f(7,:) = (ar'*sp2reg)./max(cnt'*sp2reg,1);
  for i = 1:size(sp2reg,2),
    reg = ismember(superpixels, find(sp2reg(:,i))) & ind;
    f(8,i) = median([yR(reg); 0]);
  end

  %pts = [x(superpixels(:) == 6) y(superpixels(:) == 6) z(superpixels(:) == 6)];   
  %pts = pts(~isnan(pts(:,1)),:);                                                  
  %ind = convhull(pts(:,3),pts(:,2));
  %polyarea(pts(ind,3),pts(ind,2))

  %Vertical area
  ind = angl > 90-del & angl < 90+del & ~isnan(dArea);
  ar = accumarray(superpixels(ind), dArea(ind), [nSP 1]);
  f(9,:) = ar'*sp2reg;

  %Missing area!
  missingArea = accumarray(superpixels(:), isnan(x(:)))';
  spArea = accumarray(superpixels(:), 1)';
  f(10,:) = (missingArea*sp2reg)./(spArea*sp2reg);
end


function f = areaFeaturesBatch(x, y, z, N, yR, yDir, sp2reg, superpixels)
  %Calculate the angle of each patch with the viewing direction.
  fx = 520;
  fy = 520;
  del = 15;
  nSP = max(superpixels(:));
  angl = abs(sum(N.*cat(3,x,y,z),3)./sqrt(sum(cat(3,x,y,z).^2,3)));
  %Area of the surface
  dArea = (z.*z./angl)./(fx.*fy);
  ind = find(~isnan(dArea));

  spArea = accumarray(superpixels(:), superpixels(:) > 0,[nSP 1]);
  goodSPArea = accumarray(superpixels(ind), superpixels(ind) > 0,[nSP 1]);
  f(1,:) = accumarray(superpixels(ind),dArea(ind), [nSP 1]);
  f(2,:) = f(1,:)./goodSPArea'.*spArea';
  
  %Horizontal area
  %Find the horizontal pixels and find there area
  angl = max(min(sum(bsxfun(@times, reshape(yDir,[1 1 3]), N),3),1),-1);
  angl = acosd(angl);
  
  ind = angl < 0+del & ~isnan(dArea);
  ar = accumarray(superpixels(ind),dArea(ind), [nSP 1], @sum, 0);
  f(3,:) = ar; 
  ind = angl < 0+del & ~isnan(yR);
  f(4,:) = accumarray(superpixels(ind),yR(ind),[nSP 1],@mean);
  
  ind = angl > 180-del & ~isnan(dArea);
  ar = accumarray(superpixels(ind),dArea(ind), [nSP 1], @sum, 0);
  f(5,:) = ar; 
  ind = angl > 180-del & ~isnan(yR);
  f(6,:) = accumarray(superpixels(ind),yR(ind),[nSP 1],@mean);
  
  %pts = [x(superpixels(:) == 6) y(superpixels(:) == 6) z(superpixels(:) == 6)];   
  %pts = pts(~isnan(pts(:,1)),:);                                                  
  %ind = convhull(pts(:,3),pts(:,2));
  %polyarea(pts(ind,3),pts(ind,2))

  %Vertical area
  ind = angl > 90-del & angl < 90+del & ~isnan(dArea);
  ar = accumarray(superpixels(ind), dArea(ind), [nSP 1]);
  f(7,:) = ar;
end

function f = extentFeatures(x, y, z, i, sp2reg, superpixels)
  ind = ismember(superpixels, find(sp2reg(:,i)));
  yMin = prctile(y(:),2); if(yMin > -90) yMin = -130; end
  f(1) = prctile(x(ind),5);
  f(2) = prctile(x(ind),95);
  f(3) = prctile(y(ind),5) - yMin;
  f(4) = prctile(y(ind),95) - yMin;
  f(5) = prctile(z(ind),5);
  f(6) = prctile(z(ind),95);
  f(7) = prctile(x(ind),95) - prctile(x(ind),5);
  f(8) = prctile(y(ind),95) - prctile(y(ind),5);
  f(9) = prctile(z(ind),95) - prctile(z(ind),5);
  f = f';
end
function f = planarityFeatures(pc, w, b, i, sp2reg, superpixels)
  w = w{1}; % This is the plane fit to the superpixel.
  ind = ismember(superpixels, find(sp2reg(:,i))) & ~isnan(pc(:,:,1));

  d = sum(bsxfun(@times, pc, reshape(w,[1 1 3])),3) + b;
  d = d(:);
  myStd = std(d(ind));
  myMean = mean(d(ind));
  onNum = nnz(find( d(~ind) < myMean + 3*myStd & d(~ind) > myMean - 3*myStd)) + 1;
  leftNum = nnz(find(d(~ind) < myMean - 3*myStd)) + 1;
  rightNum = nnz(find(d(~ind) > myMean + 3*myStd)) + 1;
  f = [myMean myStd log(leftNum) (onNum+leftNum)./(onNum+leftNum+rightNum) log(rightNum) (onNum+rightNum)./(onNum+leftNum+rightNum) log(leftNum./rightNum)]';

  
  ind = ismember(superpixels, find(sp2reg(:,i)));
  try
    [I J] = find(ind);
    K = convhull(I,J);
  catch
    ind = imdilate(ind,strel('disk',1));
    [I J] = find(ind);
    K = convhull(I,J);
  end
  R = roipoly(zeros(size(superpixels)),J(K),I(K)) & ~ind & ~isnan(pc(:,:,1));
  onNum = nnz(find( d(R) < myMean + 3*myStd & d(R) > myMean - 3*myStd)) + 1;
  leftNum = nnz(find(d(R) < myMean - 3*myStd)) + 1;
  rightNum = nnz(find(d(R) > myMean + 3*myStd)) + 1;
  g = [log(leftNum) (onNum+leftNum)./(onNum+leftNum+rightNum) log(rightNum) (onNum+rightNum)./(onNum+leftNum+rightNum) log(leftNum./rightNum)]';

  f = [f; g];
end
