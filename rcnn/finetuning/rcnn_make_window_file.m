function window_file = rcnn_make_window_file(imdb, out_dir, salt, image_dir, image_ext)
% rcnn_make_window_file(imdb, out_dir)
%   Makes a window file that can be used by the caffe WindowDataLayer 
%   for finetuning.
%
%   The window file format contains repeated blocks of:
%
%     # image_index 
%     img_path
%     channels 
%     height 
%     width
%     num_windows
%     class_index overlap x1 y1 x2 y2
%     <... num_windows-1 more windows follow ...>

% AUTORIGHTS
% ---------------------------------------------------------
% Copyright (c) 2014, Ross Girshick
% 
% This file is part of the R-CNN code and is available 
% under the terms of the Simplified BSD License provided in 
% LICENSE. Please retain this notice and LICENSE if you use 
% this file (or any portion of it) in your project.
% ---------------------------------------------------------

exists_or_mkdir(out_dir);
roidb = imdb.roidb_func(imdb);
window_file = fullfile(out_dir, sprintf('windowfile_%s_%s.txt', salt, imdb.name));
fid = fopen(window_file, 'wt');

channels = 3; % three channel images

for i = 1:length(imdb.image_ids)
  tic_toc_print('make window file: %d/%d\n', i, length(imdb.image_ids));
  img_path = fullfile_ext(image_dir, imdb.image_ids{i}, image_ext);
  % img_path = imdb.image_at(i);
  roi = roidb.rois(i);
  num_boxes = size(roi.boxes, 1);
  fprintf(fid, '# %d\n', i-1);
  fprintf(fid, '%s\n', img_path);
  fprintf(fid, '%d\n%d\n%d\n', ...
      channels, 0, 0);
      % imdb.sizes(i, 1), ...
      % imdb.sizes(i, 2));
  fprintf(fid, '%d\n', num_boxes);
  for j = 1:num_boxes
    [ov, label] = max(roi.overlap(j,:));
    % zero overlap => label = 0 (background)
    if ov < 1e-5
      label = 0;
      ov = 0;
    end
    bbox = roi.boxes(j,:)-1;
    fprintf(fid, '%d %.3f %d %d %d %d\n', ...
        label, ov, bbox(1), bbox(2), bbox(3), bbox(4));
  end
end

fclose(fid);
