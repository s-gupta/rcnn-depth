function [ucmThresh, meanNSP] = compute_sp_threshold(imset, numSP)
% function [ucmThresh, meanNSP] = compute_sp_threshold(imset, numSP)
  p = get_paths();
  imlist = getImageSet(imset);
  high = 1;
  low = 0.005;
  thresh = 0.005;

  while(high > low+thresh)
    mid = (low+high)/2;
    [meanNSP, nSP] = helper(imlist, mid, p);
    if(meanNSP < numSP)
      high = mid;
    else
      low = mid;
    end
    fprintf('Between %0.3f, and %0.3f\n', low, high);
  end
  
  ucmThresh = mid;
end


function [meanNSP, nSP] = helper(imlist, ucm_thresh, p)
  parfor i = 1:length(imlist),
    imname = imlist{i};
    dt = load(fullfile_ext(p.ucm_dir, 'multi', imname, 'mat'));ucm2 = dt.ucm2;
    sp = bwlabel(dt.ucm2 < ucm_thresh); sp = sp(2:2:end, 2:2:end);
    nSP(i) = max(sp(:));
  end
  meanNSP = mean(nSP);
end
