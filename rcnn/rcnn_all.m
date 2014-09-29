function res = rcnn_all(task, model_typ, trainset, testset)
  global RUNNAME
  RUNNAME = 'release';
  p = get_paths(RUNNAME);
  c = benchmarkPaths();
  NYU_ROOT_DIR = c.dataDir;
  REGIONDIR = fullfile(p.box_dir, 'release');
  SALT = 'release';
  MAX_BOXES = 2000;

  imdb = imdb_from_nyud2(NYU_ROOT_DIR, trainset, task, REGIONDIR, SALT, MAX_BOXES);
  imdb.roidb_func = @roidb_from_nyud2;
  
  feat_opts = struct('featDir', '');
  feat_opts(1) = [];
  global RCNN_CONFIG_OVERRIDE;
  switch model_typ
    case 'rgb_hha',
      conf_override.sub_dir = sprintf('rgb_hha_%d_%s', 30000, trainset);
      feat_opts(end+1) = struct('featDir', fullfile(p.cnnF_cache_dir, 'rgb_30000', imdb.dataset_name));
      feat_opts(end+1) = struct('featDir', fullfile(p.cnnF_cache_dir, 'hha_30000', imdb.dataset_name));
    case 'rgb',
      conf_override.sub_dir = sprintf('rgb_%d_%s', 30000, trainset);
      feat_opts(end+1) = struct('featDir', fullfile(p.cnnF_cache_dir, 'rgb_30000', imdb.dataset_name));
    case 'hha',
      conf_override.sub_dir = sprintf('hha_%d_%s', 30000, trainset);
      feat_opts(end+1) = struct('featDir', fullfile(p.cnnF_cache_dir, 'hha_30000', imdb.dataset_name));
    case 'disparity',
      conf_override.sub_dir = sprintf('disparity_%d_%s', 30000, trainset);
      feat_opts(end+1) = struct('featDir', fullfile(p.cnnF_cache_dir, 'disparity_30000', imdb.dataset_name));
  end
  RCNN_CONFIG_OVERRIDE = @() conf_override;
  rcnn_train_opts = {'feat_opts', feat_opts};
  
  rcnn_models = rcnn_train(imdb, rcnn_train_opts{:});
  
  if(isstr(testset)), testset = {testset}; end
  for i = 1:length(testset),
    imdb = imdb_from_nyud2(NYU_ROOT_DIR, testset{i}, task, REGIONDIR, SALT, MAX_BOXES);
    imdb.roidb_func = @roidb_from_nyud2;
    res{i} = rcnn_test(rcnn_models, imdb, '');
  end
end
