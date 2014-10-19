function out = getProbCube(outp, imName, numClass)
% function out = getProbCube(outp, i, numClass)
  dt = load(fullfile(outp, strcat(imName, '.mat')));
  for j = 1:size(dt.scores, 1), 
    si = dt.scores(j,:); 
    out(:,:,j) = si(dt.superpixels); 
  end
end
