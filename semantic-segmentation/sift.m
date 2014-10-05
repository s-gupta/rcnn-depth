function bow_map = sift(I, siftParam, vocab)
  paths = getPaths();
  
  % Sift Param
  % siftParam = struct('ds_sampling', 1, 'scales', 1.2, 'descriptor', 'opponentsift');
  
  %% Compute the SIFT Features
  % Generate three filenames
  ext = {'png', 'txt'};
  r = randsample(100000, length(ext), false);
  pid = getPID();
  for i = 1:length(ext),
    %f{i} = fullfile('/tmp', sprintf('sgupta-sift-%07d-%06d.%s', pid, r(i)), ext{i});
    f{i} = fullfile('/dev/shm', sprintf('sgupta-sift-%07d-%06d.%s', pid, r(i), ext{i}));
  end

  imwrite(im2uint8(I), f{1}, 'png');
  str = sprintf('%s %s --detector densesampling --ds_spacing %d --ds_scales %0.4f --descriptor %s --output %s\n', paths.siftLib, f{1}, siftParam.ds_sampling, siftParam.scales, siftParam.descriptor, f{2});
  system(str);
  
  info = imfinfo(f{1});
  feats = feat_txt2mat(f{2}, siftParam.scales);

  % Remove the temporary files.
  for i = 1:length(f)
    system(sprintf('rm %s &', f{i}));
  end


  %% Compute the SIFT bow here
  data_img = feats.data';
  feats.pos = double(feats.pos); 
  pos = sub2ind([info.Height, info.Width], feats.pos(:,1), feats.pos(:,2));
  cnt_id = 1:numel(pos);
  im_pos = zeros([info.Height info.Width]);
  im_pos(pos) = cnt_id;
  clear feats;
  
  [drop, binsa] = min(vl_alldist(vocab, single(data_img)), [], 1) ;
  
  ind = (im_pos > 0);
  
  bow_map = zeros(size(im_pos));
  bow_map(ind) = binsa(im_pos(ind));
end
