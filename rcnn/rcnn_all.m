function res = rcnn_all(task, model_typ, region_features, region_task, trainset, testset, imdb_salt, salt, region_training, varargin)
% function res = rcnn_all(task, model_typ, region_features, region_task, trainset, testset, imdb_salt, salt, region_training, varargin)
  % global RUNNAME
  % RUNNAME = 'release';
  % p = get_paths(RUNNAME);

  p = get_paths();
  c = benchmarkPaths();
  NYU_ROOT_DIR = c.dataDir;
  REGIONDIR = fullfile(p.output_dir, 'regions', 'release-gt-inst');
  SALT = imdb_salt;
  MAX_BOXES = 2000;

  imdb = imdb_from_nyud2(NYU_ROOT_DIR, trainset, task, REGIONDIR, SALT, MAX_BOXES);
  if region_task
    imdb.roidb_func = @roidb_from_nyud2_region;
    imdb.eval_func = @imdb_eval_nyud2_region;
    box_or_region = 'region';
  else
    imdb.eval_func = @imdb_eval_nyud2;
    imdb.roidb_func = @roidb_from_nyud2;
    box_or_region = 'box';
  end

  if region_task, assert(region_features == true, 'region_fetures must be true if region_task is true'); end
  
  feat_opts = struct('featDir', '');
  feat_opts(1) = [];
  global RCNN_CONFIG_OVERRIDE;
  switch model_typ
    case 'rgb_hha',
      conf_override.sub_dir = sprintf('rgb_hha_%d_%s', 30000, trainset);
      feat_opts(end+1) = struct('featDir', fullfile(p.cnnF_cache_dir, 'rgb_30000', imdb.dataset_name));
      feat_opts(end+1) = struct('featDir', fullfile(p.cnnF_cache_dir, 'hha_30000', imdb.dataset_name));
      if(region_features)
        feat_opts(end+1) = struct('featDir', fullfile(p.cnnF_cache_dir, 'rgb_region_30000', imdb.dataset_name));
        feat_opts(end+1) = struct('featDir', fullfile(p.cnnF_cache_dir, 'hha_region_30000', imdb.dataset_name));
      end
    case 'rgb',
      conf_override.sub_dir = sprintf('rgb_%d_%s', 30000, trainset);
      feat_opts(end+1) = struct('featDir', fullfile(p.cnnF_cache_dir, 'rgb_30000', imdb.dataset_name));
      if region_features
        feat_opts(end+1) = struct('featDir', fullfile(p.cnnF_cache_dir, 'rgb_region_30000', imdb.dataset_name));
      end
    case 'hha',
      conf_override.sub_dir = sprintf('hha_%d_%s', 30000, trainset);
      feat_opts(end+1) = struct('featDir', fullfile(p.cnnF_cache_dir, 'hha_30000', imdb.dataset_name));
      if region_features
        feat_opts(end+1) = struct('featDir', fullfile(p.cnnF_cache_dir, 'hha_region_30000', imdb.dataset_name));
      end
    case 'disparity',
      conf_override.sub_dir = sprintf('disparity_%d_%s', 30000, trainset);
      feat_opts(end+1) = struct('featDir', fullfile(p.cnnF_cache_dir, 'disparity_30000', imdb.dataset_name));
      if region_features,
        feat_opts(end+1) = struct('featDir', fullfile(p.cnnF_cache_dir, 'disparity_region_30000', imdb.dataset_name));
      end
  end
  conf_override.sub_dir = sprintf('%s_region-features-%d_region-task-%d%s', ...
    conf_override.sub_dir, region_features, region_task, salt);

  if(region_task),
    n = length(feat_opts)/2; feat_opts = feat_opts([(n+1):end, 1:n]); end

  RCNN_CONFIG_OVERRIDE = @() conf_override;
  rcnn_train_opts = {'feat_opts', feat_opts};
  if region_task && region_training,
    rcnn_models = rcnn_region_train(imdb, rcnn_train_opts{:}, varargin{:});
  else
    rcnn_models = rcnn_train(imdb, rcnn_train_opts{:}, varargin{:});
  end

  
  if(isstr(testset)), testset = {testset}; end
  for i = 1:length(testset),
    imdb = imdb_from_nyud2(NYU_ROOT_DIR, testset{i}, task, REGIONDIR, SALT, MAX_BOXES);
    if region_task
      imdb.roidb_func = @roidb_from_nyud2_region;
      imdb.eval_func = @imdb_eval_nyud2_region;
      box_or_region = 'region';
    else
      imdb.eval_func = @imdb_eval_nyud2;
      imdb.roidb_func = @roidb_from_nyud2;
      box_or_region = 'box';
    end
    res{i} = rcnn_test(rcnn_models, imdb, '', box_or_region);
  end
end
