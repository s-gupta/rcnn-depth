if strcmp(jobName, 'save_disparity')
  p = get_paths();
  C = cropCamera(getCameraParam('color'));
  outDir = p.ft_disparity_dir;
  imList = getImageSet('all');
  parfor i = 1:length(imList),
    args{i} = {imList{i}, C, outDir, [], []}; 
    saveDisparity(args{i}{:});
  end
end

if strcmp(jobName, 'save_color')
  p = get_paths();
  outDir = p.ft_image_dir;
  imList = getImageSet('all');
  parfor i = 1:length(imList),
    I = getImage(imList{i}, 'images');
    imwrite(I, fullfile_ext(outDir, imList{i}, 'png'), 'png');
  end
end

if strcmp(jobName, 'save_hha')
  p = get_paths();
  C = cropCamera(getCameraParam('color'));
  outDir = p.ft_hha_dir;
  imList = getImageSet('all');
  parfor i = 1:length(imList),
    args{i} = {imList{i}, C, outDir, [], []}; 
    saveHHA(args{i}{:});
  end
end

if strcmp(jobName, 'write_window_file')
  p = get_paths();
  c = benchmarkPaths();
  NYU_ROOT_DIR = c.dataDir;
  REGIONDIR = fullfile(p.output_dir, 'regions', 'release');
  SALT = 'release';
  MAX_BOXES = 2000;
  task = 'task-detection-with-cabinet';
  
  trainset = 'train';
  imdb = imdb_from_nyud2(c.dataDir, trainset, task, REGIONDIR, SALT, MAX_BOXES);
  imdb.roidb_func = @roidb_from_nyud2;
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



if strcmp(jobName, 'hha_cache_features')
  p = get_paths();
  c = benchmarkPaths();
  NYU_ROOT_DIR = c.dataDir;
  REGIONDIR = fullfile(p.output_dir, 'regions', 'release');
  SALT = 'release';
  MAX_BOXES = 2000;
  task = 'task-detection-with-cabinet';

  imset = 'all';
  imdb = imdb_from_nyud2(NYU_ROOT_DIR, imset, task, REGIONDIR, SALT, MAX_BOXES);
  imdb.roidb_func = @roidb_from_nyud2;
 
    
  p = get_paths();
  image_dir = fullfile(p.ft_hha_dir);
  image_ext = 'png';
  snapshot = 30000;
  net_file = fullfile_ext(p.snapshot_dir, sprintf('nyud2_finetune_hha_iter_%d', snapshot), 'caffemodel');
  feat_cache_dir = p.cnnF_cache_dir;
  net_def_file = fullfile('nyud2_finetuning', 'imagenet_hha_256_fc6.prototxt');
  mean_file = fullfile_ext(p.mean_file_hha, 'mat');
  cache_name = sprintf('hha_%d', snapshot);
  args = {};
  st = 1; sp = 1; e = 0; gpu_id = 1;
  args{1} = {'start', st, 'step', sp, 'end', e, ...
    'image_dir', image_dir, 'image_ext', image_ext, ...
    'feat_cache_dir', feat_cache_dir, ...
    'net_def_file', net_def_file, 'net_file', net_file, 'mean_file', mean_file, ...
    'cache_name', cache_name, 'gpu_id', gpu_id};
  rcnn_cache_features(imdb, args{1}{:});
end

if strcmp(jobName, 'disparity_cache_features')
  p = get_paths();
  c = benchmarkPaths();
  NYU_ROOT_DIR = c.dataDir;
  REGIONDIR = fullfile(p.output_dir, 'regions', 'release');
  SALT = 'release';
  MAX_BOXES = 2000;
  task = 'task-detection-with-cabinet';

  imset = 'all';
  imdb = imdb_from_nyud2(NYU_ROOT_DIR, imset, task, REGIONDIR, SALT, MAX_BOXES);
  imdb.roidb_func = @roidb_from_nyud2;
  
  p = get_paths();
  image_dir = fullfile(p.ft_disparity_dir);
  image_ext = 'png';
  snapshot = 30000;
  net_file = fullfile_ext(p.snapshot_dir, sprintf('nyud2_finetune_disparity_iter_%d', snapshot), 'caffemodel');
  feat_cache_dir = p.cnnF_cache_dir;
  net_def_file = fullfile('nyud2_finetuning', 'imagenet_disparity_256_fc6.prototxt');
  mean_file = fullfile_ext(p.mean_file_disparity, 'mat');
  cache_name = sprintf('disparity_%d', snapshot);
  args = {};
  st = 1; sp = 1; e = 0; gpu_id = 0;
  args{1} = {'start', st, 'step', sp, 'end', e, ...
    'image_dir', image_dir, 'image_ext', image_ext, ...
    'feat_cache_dir', feat_cache_dir, ...
    'net_def_file', net_def_file, 'net_file', net_file, 'mean_file', mean_file, ...
    'cache_name', cache_name, 'gpu_id', gpu_id};
  rcnn_cache_features(imdb, args{1}{:});
end



if strcmp(jobName, 'color_cache_features')
  p = get_paths();
  c = benchmarkPaths();
  NYU_ROOT_DIR = c.dataDir;
  REGIONDIR = fullfile(p.output_dir, 'regions', 'release');
  SALT = 'release';
  MAX_BOXES = 2000;
  task = 'task-detection-with-cabinet';

  imset = 'all';
  imdb = imdb_from_nyud2(NYU_ROOT_DIR, imset, task, REGIONDIR, SALT, MAX_BOXES);
  imdb.roidb_func = @roidb_from_nyud2;
  
  p = get_paths();
  image_dir = fullfile(p.ft_image_dir);
  image_ext = 'png';
  snapshot = 30000;
  net_file = fullfile_ext(p.snapshot_dir, sprintf('nyud2_finetune_color_iter_%d', snapshot), 'caffemodel');
  feat_cache_dir = p.cnnF_cache_dir;
  net_def_file = fullfile('nyud2_finetuning', 'imagenet_color_256_fc6.prototxt');
  mean_file = fullfile_ext(p.mean_file_color, 'mat');
  cache_name = sprintf('rgb_%d', snapshot);
  args = {};
  st = 1; sp = 1; e = 0; gpu_id = 0;
  args{1} = {'start', st, 'step', sp, 'end', e, ...
    'image_dir', image_dir, 'image_ext', image_ext, ...
    'feat_cache_dir', feat_cache_dir, ...
    'net_def_file', net_def_file, 'net_file', net_file, 'mean_file', mean_file, ...
    'cache_name', cache_name, 'gpu_id', gpu_id};
  rcnn_cache_features(imdb, args{1}{:});
end

if strcmp(jobName, 'rcnn_train'),
  % The final model
  res = rcnn_all('task-detection', 'rgb_hha', 'trainval', 'test');

  % Control experiments
  % res = rcnn_all('task-detection', 'rgb_hha', 'train', 'val');
  % res = rcnn_all('task-detection', 'hha', 'train', 'val');
  % res = rcnn_all('task-detection', 'rgb', 'train', 'val');
  % res = rcnn_all('task-detection', 'disparity', 'train', 'val');
end

if strcmp(jobName, 'proto_to_mat_mean')
  p = get_paths();
  image_mean = caffe('read_mean', fullfile_ext(p.mean_file_color, 'proto'));
  image_mean = permute(image_mean, [2 1 3]);
  save(fullfile_ext(p.mean_file_color, 'mat'), 'image_mean');
  
  image_mean = caffe('read_mean', fullfile_ext(p.mean_file_hha, 'proto'));
  image_mean = permute(image_mean, [2 1 3]);
  save(fullfile_ext(p.mean_file_hha, 'mat'), 'image_mean');
  
  image_mean = caffe('read_mean', fullfile_ext(p.mean_file_disparity, 'proto'));
  image_mean = permute(image_mean, [2 1 3]);
  save(fullfile_ext(p.mean_file_disparity, 'mat'), 'image_mean');
end
