args = {};

if strcmp(jobName, 'region_write_window_file')
  p = get_paths(); c = benchmarkPaths();
  NYU_ROOT_DIR = c.dataDir; REGIONDIR = fullfile(p.output_dir, 'regions', 'release-gt-inst');
  SALT = 'release-regions';
  MAX_BOXES = 2000;
  task = 'task-detection-with-cabinet';
  
  trainset = 'test';
  imdb = imdb_from_nyud2(c.dataDir, trainset, task, REGIONDIR, SALT, MAX_BOXES);
  imdb.roidb_func = @roidb_from_nyud2_region;
  roidb = imdb.roidb_func(imdb);
  spdir = fullfile(p.detection_dir,  'sp');
  sp2regdir = fullfile(p.detection_dir,  'sp2reg')
  write_superpixels_sp2reg(imdb, roidb, spdir, sp2regdir);
  
  imdirs = {}; imexts = {'png', 'png', 'png'};
  imlist = imdb.image_ids; channels = 3;
  list = {p.ft_hha_dir, sp2regdir, spdir};
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

  write_window_file(imdirs, imexts, imlist, channels, list, window_file)
  
  window_file_train_color = rcnn_make_window_file(imdb, fullfile(p.wf_dir), 'color', p.ft_image_dir, 'png');
  window_file_train_hha = rcnn_make_window_file(imdb, fullfile(p.wf_dir), 'hha', p.ft_hha_dir, 'png');

  keyboard;

  trainset = 'val';
  imdb = imdb_from_nyud2(NYU_ROOT_DIR, trainset, task, REGIONDIR, SALT, MAX_BOXES);
  imdb.roidb_func = @roidb_from_nyud2;
  window_file_val_color = rcnn_make_window_file(imdb, fullfile(p.wf_dir), 'color', p.ft_image_dir, 'png');
  window_file_val_hha = rcnn_make_window_file(imdb, fullfile(p.wf_dir), 'hha', p.ft_hha_dir, 'png');
  window_file_val_disparity = rcnn_make_window_file(imdb, fullfile(p.wf_dir), 'disparity', p.ft_disparity_dir, 'png');
  
  file_name = {'nyud2_finetune_solver.color', 'nyud2_finetune_solver.hha', 'nyud2_finetune_solver.disparity', ...
  'nyud2_finetune_train.color', 'nyud2_finetune_train.hha', 'nyud2_finetune_train.disparity', ...
  'nyud2_finetune_val.color', 'nyud2_finetune_val.hha', 'nyud2_finetune_val.disparity', ...
  'imagenet_deploy.prototxt'};

  args = {'PROTODIR', p.proto_dir, ...
    'WINDOW_FILE_TRAIN_COLOR', window_file_train_color, ...
    'WINDOW_FILE_TRAIN_HHA', window_file_train_hha, ...
    'WINDOW_FILE_TRAIN_DISPARITY', window_file_train_disparity, ...
    'WINDOW_FILE_VAL_COLOR', window_file_val_color, ... 
    'WINDOW_FILE_VAL_HHA', window_file_val_hha, ...
    'WINDOW_FILE_VAL_DISPARITY', window_file_val_disparity, ...
    'MEAN_FILE_COLOR', fullfile_ext(p.mean_file_color, 'proto'), ...
    'MEAN_FILE_HHA', fullfile_ext(p.mean_file_hha, 'proto'), ...
    'MEAN_FILE_DISPARITY', fullfile_ext(p.mean_file_disparity, 'proto'), ...
    'N_CLASSES', sprintf('%d', imdb.num_classes+1), ...
    'SNAPSHOTDIR', p.snapshot_dir,...
  };
  make_protofiles(p.proto_dir, file_name, args);

  % Generate the caffe command
  fprintf('GLOG_logstdtoerr=1 caffe/build/tools/caffe.bin train -gpu 0 -iterations 30000 -solver %s -weights %s -model %s/imagenet_deploy.prototxt 2>&1 | tee %s\n\n\n', fullfile(p.proto_dir, 'nyud2_finetune_solver.color'), p.caffe_net, p.proto_dir, fullfile(p.proto_dir, 'color.log'));
  fprintf('GLOG_logstdtoerr=1 caffe/build/tools/caffe.bin train -gpu 0 -iterations 30000 -solver %s -weights %s -model %s/imagenet_deploy.prototxt 2>&1 | tee %s\n\n\n', fullfile(p.proto_dir, 'nyud2_finetune_solver.hha'), p.caffe_net, p.proto_dir, fullfile(p.proto_dir, 'hha.log'));

  % Control experiment for finetuning on the disparity only..
  % fprintf('GLOG_logstdtoerr=1 caffe/build/tools/caffe.bin train -gpu 1 -iterations 30000 -solver %s -weights %s -model %s/imagenet_deploy.prototxt 2>&1 | tee %s\n\n\n', fullfile(p.proto_dir, 'nyud2_finetune_solver.disparity'), p.caffe_net, p.proto_dir, fullfile(p.proto_dir, 'disparity.log'));
end

if strcmp(jobName, 'vis_regions')
  p = get_paths(); c = benchmarkPaths();
  NYU_ROOT_DIR = c.dataDir; REGIONDIR = fullfile(p.output_dir, 'regions', 'release-gt-inst');
  SALT = 'release-regions';
  MAX_BOXES = 2000;
  task = 'task-detection-with-cabinet';

  trainset = 'img_6449';
  imdb = imdb_from_nyud2(c.dataDir, trainset, task, REGIONDIR, SALT, MAX_BOXES);
  imdb.roidb_func = @roidb_from_nyud2_region;

  roidb = imdb.roidb_func(imdb);
  roi = roidb.rois;
  clss = 5; %[1:20];
  I = getImage(imdb.image_ids{1}, 'images');
  figure(1); imagesc(I);
  for i = clss,
    figure(2);
    [~, ind] = sort(roi.overlap(:,i), 'descend');
    for j = 1:1:100,
      sp2regi = roi.sp2reg(ind(j),:);
      imagesc(sp2regi(roi.sp));
      title(sprintf('%s %0.3f', imdb.classes{i}, roi.overlap(ind(j),i)));
      pause;
    end
  end
end

