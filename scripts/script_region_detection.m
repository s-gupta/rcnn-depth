args = {};

if strcmp(jobName, 'region_write_window_file')
  p = get_paths(); c = benchmarkPaths();
  NYU_ROOT_DIR = c.dataDir; REGIONDIR = fullfile(p.output_dir, 'regions', 'release-gt-inst');
  SALT = 'release-regions';
  MAX_BOXES = 2000;
  task = 'task-detection-with-cabinet';
  
  trainset = 'train';
  imdb = imdb_from_nyud2(c.dataDir, trainset, task, REGIONDIR, SALT, MAX_BOXES);
  imdb.roidb_func = @roidb_from_nyud2_region;


  window_file_train_color = rcnn_make_window_file(imdb, fullfile(p.wf_dir), 'color', p.ft_image_dir, 'png');
  window_file_train_hha = rcnn_make_window_file(imdb, fullfile(p.wf_dir), 'hha', p.ft_hha_dir, 'png');
  window_file_train_disparity = rcnn_make_window_file(imdb, fullfile(p.wf_dir), 'disparity', p.ft_disparity_dir, 'png');

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



