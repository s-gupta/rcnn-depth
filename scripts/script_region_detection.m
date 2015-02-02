args = {};

if strcmp(jobName, 'region_write_window_file')
  p = get_paths(); c = benchmarkPaths();
  NYU_ROOT_DIR = c.dataDir; REGIONDIR = fullfile(p.output_dir, 'regions', 'release-gt-inst');
  SALT = 'region';
  MAX_BOXES = 2000;
  task = 'task-detection-with-cabinet';
  
  spdir = fullfile(p.detection_dir,  'sp');
  sp2regdir = fullfile(p.detection_dir,  'sp2reg');
  imexts = {'png', 'png', 'png'};
  imsets = {'trainval', 'test'}; 

  channels = 3;
  for i = 1:2,
    imset = imsets{i}; 
    imdb = imdb_from_nyud2(c.dataDir, imset, task, REGIONDIR, SALT, MAX_BOXES);
    % imdb.roidb_func = @roidb_from_nyud2_region;
    imdb.roidb_func = @roidb_from_nyud2;
    roidb = imdb.roidb_func(imdb);
    % write_superpixels_sp2reg(imdb, roidb, spdir, sp2regdir);
    list = {}; 
    for i = 1:length(roidb.rois),
      roi = roidb.rois(i); [ov, label] = max(roi.overlap, [], 2);
      ov(ov < 1e-5) = 0; label(ov < 1e-5) = 0; list{i} = cat(2, label, ov, roi.boxes-1);
    end
    
    imlist = imdb.image_ids;
    imdirs = {p.ft_hha_dir, spdir, sp2regdir};
    window_file = fullfile(p.detection_dir, 'finetuning', 'v3', 'wf', sprintf('ft_box_hha_%s', imset)); 
    write_window_file(imdirs, imexts, imlist, channels, list, window_file);

    imdirs = {p.ft_image_dir, spdir, sp2regdir};
    window_file = fullfile(p.detection_dir, 'finetuning', 'v3', 'wf', sprintf('ft_box_rgb_%s', imset)); 
    write_window_file(imdirs, imexts, imlist, channels, list, window_file);
  end
end

%   GLOG_logtostderr=1 ../../caffe/build_region_vader/tools/caffe.bin train \
%   -gpu 1 \
%   -model ../eccv14-cachedir/release/detection/finetuning/v3/proto/finetune_box.hha \
%   -weights caffe-data/caffe_reference_imagenet_model \
%   -solver ../eccv14-cachedir/release/detection/finetuning/v3/proto/solver_box.hha 2>&1 \
%   | tee ../eccv14-cachedir/release/detection/finetuning/v3/log_box.hha
% 
%   GLOG_logtostderr=1 ../../caffe/build_region_vader/tools/caffe.bin train \
%   -gpu 0 \
%   -model ../eccv14-cachedir/release/detection/finetuning/v3/proto/finetune_box.rgb \
%   -weights caffe-data/caffe_reference_imagenet_model \
%   -solver ../eccv14-cachedir/release/detection/finetuning/v3/proto/solver_box.rgb 2>&1 \
%   | tee ../eccv14-cachedir/release/detection/finetuning/v3/log_box.rgb
%  
%   GLOG_logtostderr=1 ../../caffe/build_region_vader/tools/extract_features \
%   ../eccv14-cachedir/release/detection/finetuning/v3/snapshot/rgb_iter_30000.caffemodel \
%   ../eccv14-cachedir/release/detection/finetuning/v3/proto/finetune_region.rgb \
%   hdf5 fc6 \
%   ../eccv14-cachedir/release/detection/finetuning/v3/features/fc6_30000_region_trainval.rgb.h5 \
%   13000 GPU 0 \
%   | tee ../eccv14-cachedir/release/detection/finetuning/v3/fe_log_region_trainval.rgb
%
%   GLOG_logtostderr=1 ../../caffe/build_region_vader/tools/extract_features \
%   ../eccv14-cachedir/release/detection/finetuning/v3/snapshot/hha_iter_30000.caffemodel \
%   ../eccv14-cachedir/release/detection/finetuning/v3/proto/finetune_region.hha \
%   hdf5 fc6 \
%   ../eccv14-cachedir/release/detection/finetuning/v3/features/fc6_30000_region_trainval.hha.h5 \
%   13000 GPU 0 \
%   | tee ../eccv14-cachedir/release/detection/finetuning/v3/fe_log_region_trainval.hha




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

if strcmp(jobName, 'check_boxes'),
  imset = 'train';
  imlist = getImageSet(imset);
  pt = load(['cache/release/detection/finetuning/v1/wf/ft_hha_' imset '.mat']);
  for i = 1:length(imlist),
    dt = load(fullfile_ext('cache/release/detection/feat_cache/hha_30000/nyud2_release/', imlist{i}, 'mat'), 'boxes');
    assert(isequal(dt.boxes, single(pt.list{i}(:, 3:6))+1));
  end
end

if strcmp(jobName, 'hha_cache_region_features')
  p = get_paths(); c = benchmarkPaths();
  NYU_ROOT_DIR = c.dataDir;
  REGIONDIR = fullfile(p.output_dir, 'regions', 'release-gt-inst');
  SALT = 'release'; MAX_BOXES = 2000;
  task = 'task-detection-with-cabinet';

  imset = 'test';
  imdb = imdb_from_nyud2(NYU_ROOT_DIR, imset, task, REGIONDIR, SALT, MAX_BOXES);
  imdb.roidb_func = @roidb_from_nyud2_region;
    
  image_dir = fullfile(p.ft_hha_dir); hha_or_rgb = 'hha'; 
  
  image_ext = 'png'; snapshot = 30000;
  net_file = fullfile_ext(p.snapshot_dir, sprintf('nyud2_finetune_region_%s_iter_%d', hha_or_rgb, snapshot), 'caffemodel');
  feat_cache_dir = p.cnnF_cache_dir;
  net_def_file = fullfile('nyud2_finetuning', 'imagenet_hha_256_fc6.prototxt');
  mean_file = fullfile_ext(p.mean_file_hha, 'mat');
  cache_name = sprintf('hha_region_%d', snapshot);
  args = {};
  % st = 1; sp = 1; e = 0; gpu_id = 1;
  args{1} = {'start', st, 'step', sp, 'end', e, ...
    'image_dir', image_dir, 'image_ext', image_ext, ...
    'feat_cache_dir', feat_cache_dir, ...
    'net_def_file', net_def_file, 'net_file', net_file, 'mean_file', mean_file, ...
    'cache_name', cache_name, 'gpu_id', gpu_id};
  rcnn_cache_features(imdb, true, args{1}{:});
end

if strcmp(jobName, 'rgb_cache_region_features')
  p = get_paths(); c = benchmarkPaths();
  NYU_ROOT_DIR = c.dataDir;
  REGIONDIR = fullfile(p.output_dir, 'regions', 'release-gt-inst');
  SALT = 'release'; MAX_BOXES = 2000;
  task = 'task-detection-with-cabinet';

  % imset = 'train';
  imdb = imdb_from_nyud2(NYU_ROOT_DIR, imset, task, REGIONDIR, SALT, MAX_BOXES);
  imdb.roidb_func = @roidb_from_nyud2_region;
    
  image_dir = fullfile(p.ft_image_dir); hha_or_rgb = 'rgb'; 
  
  image_ext = 'png'; snapshot = 30000;
  net_file = fullfile_ext(p.snapshot_dir, sprintf('nyud2_finetune_region_%s_iter_%d', hha_or_rgb, snapshot), 'caffemodel');
  feat_cache_dir = p.cnnF_cache_dir;
  net_def_file = fullfile('nyud2_finetuning', 'imagenet_color_256_fc6.prototxt');
  mean_file = fullfile_ext(p.mean_file_color, 'mat');
  cache_name = sprintf('rgb_region_%d', snapshot);
  args = {};
  % st = 1; sp = 1; e = 0; gpu_id = 1;
  args{1} = {'start', st, 'step', sp, 'end', e, ...
    'image_dir', image_dir, 'image_ext', image_ext, ...
    'feat_cache_dir', feat_cache_dir, ...
    'net_def_file', net_def_file, 'net_file', net_file, 'mean_file', mean_file, ...
    'cache_name', cache_name, 'gpu_id', gpu_id};
  rcnn_cache_features(imdb, true, args{1}{:});
end

if strcmp(jobName, 'fe'),
  p = get_paths(); c = benchmarkPaths();
  NYU_ROOT_DIR = c.dataDir;
  REGIONDIR = fullfile(p.output_dir, 'regions', 'release-gt-inst');
  SALT = 'release'; MAX_BOXES = 2000;
  task = 'task-detection-with-cabinet';
  imsets = {'trainval', 'test'};
  box_or_region = {'region', 'box'};
  type = {'rgb', 'hha'};
  for i = 2,
    for j = 1:2,
      for k = [2, 1],
        imset = imsets{i};
        imdb = imdb_from_nyud2(NYU_ROOT_DIR, imset, task, REGIONDIR, SALT, MAX_BOXES);
        if(strcmp(box_or_region{j}, 'region'))
          imdb.roidb_func = @roidb_from_nyud2_region;
        elseif(strcmp(box_or_region{j}, 'box'))
          imdb.roidb_func = @roidb_from_nyud2;
        end

        % h5_file = sprintf('cache/release/detection/feat_cache/rgb_region_30000/%s.h5', imset);
        h5_file = sprintf('/data1/sgupta/tmptmp/nyud2/fc6_30000_%s_%s.%s.h5', box_or_region{j}, imset, type{k})
        window_file = sprintf('cache/release/detection/finetuning/v3/wf/ft_%s_%s_%s.mat', ...
          box_or_region{j}, type{k}, imset);
        output_dir = sprintf('cache/release_ft-trainval/detection/feat_cache/%s_%s_30000/nyud2_release/', ...
          type{k}, box_or_region{j})
        mkdir_if_missing(output_dir);
        h5_to_mat_fast(h5_file, imdb, window_file, output_dir);
      end
    end
  end 
end

if strcmp(jobName, 'hha_cache_features')
  global RUNNAME
  RUNNAME = 'release_ft-trainval';
  p = get_paths();
  c = benchmarkPaths();
  NYU_ROOT_DIR = c.dataDir;
  REGIONDIR = fullfile(p.output_dir, 'regions', 'release-gt-inst');
  SALT = 'matlab';
  MAX_BOXES = 2000;
  task = 'task-detection-with-cabinet';

  % imset = 'trainval'; type = 'hha'; gpu_id = 0;
  imdb = imdb_from_nyud2(NYU_ROOT_DIR, imset, task, REGIONDIR, SALT, MAX_BOXES);
  
  % box_or_region = 'box'; imdb.roidb_func = @roidb_from_nyud2; region = false;
  box_or_region = 'region'; imdb.roidb_func = @roidb_from_nyud2_region; region = true;
 
    
  p = get_paths();
  image_dir = fullfile(p.(sprintf('ft_%s_dir', type)));
  image_ext = 'png';
  snapshot = 30000;
  net_file = fullfile_ext(p.snapshot_dir, sprintf('%s_%s_iter_%d', box_or_region, type, snapshot), 'caffemodel');
  feat_cache_dir = p.cnnF_cache_dir;
  net_def_file = fullfile('nyud2_finetuning', sprintf('imagenet_%s_256_fc6.prototxt', type));
  mean_file = fullfile_ext(p.(sprintf('mean_file_%s', type)), 'mat');
  cache_name = sprintf('%s_%s_%d', type, box_or_region, snapshot);
  args = {};
  st = 1; sp = 1; e = 0;
  args{1} = {'start', st, 'step', sp, 'end', e, ...
    'image_dir', image_dir, 'image_ext', image_ext, ...
    'feat_cache_dir', feat_cache_dir, ...
    'net_def_file', net_def_file, 'net_file', net_file, 'mean_file', mean_file, ...
    'cache_name', cache_name, 'gpu_id', gpu_id};
  rcnn_cache_features(imdb, region, args{1}{:});
end

% jobName = 'hha_cache_features'; type = 'hha'; gpu_id = 0; imset = 'trainval'; script_region_detection; exit; 
% jobName = 'hha_cache_features'; type = 'hha'; gpu_id = 0; imset = 'test'; script_region_detection; exit;
% jobName = 'hha_cache_features'; type = 'rgb'; gpu_id = 1; imset = 'trainval'; script_region_detection; exit; 
% jobName = 'hha_cache_features'; type = 'rgb'; gpu_id = 1; imset = 'test'; script_region_detection; exit;



% Check the extracted data..
if strcmp(jobName, 'feat_norm_diff')
  dbstop in rcnn_features at 32
  imset = 'img_5003';
  st = 1; sp = 1; e = 0; gpu_id = 0; jobName = 'rgb_cache_region_features'; script_region_detection
  a = read_h5_file('cache/release/detection/feat_cache/rgb_region_30000/data.h5', 'data-[0-9]*');

  batches{1} = batches{1}(:,:,:,1:128);
  for i = 1:128, 
    im1 = uint8(permute(a{1}(:,:,1:3,i)+128, [2 1 3]));
    im2 = uint8(permute(batches{1}(:,:,1:3,i)+128, [2 1 3]));
    figure(1); 
    subplot(1,3,1); imagesc(im1); subplot(1,3,2); imagesc(im2);
    subplot(1,3,3); imagesc(sqrt(sum((im1-im2).^2, 3))); colormap jet; colorbar; title(norm(double(im1(:))-double(im2(:)))./norm(double(im1(:))));
    pause; 
  end

  b1 = cat(4, batches{1}, batches{1}); f1 = caffe('forward', {b1}); 
  b2 = cat(4, a{1}, a{1}); f2 = caffe('forward', {b2}); 
  f1 = f1{1}(:); f2 = f2{1}(:);
  norm(f1-f2)./norm(f2);
end

if strcmp(jobName, 'rcnn')
  res_tr{1} = rcnn_all('task-detection', 'rgb', 0, 0, 'train', 'val'); 
  res_tr{end+1} = rcnn_all('task-detection', 'hha', 0, 0, 'train', 'val'); 
  res_tr{end+1} = rcnn_all('task-detection', 'rgb_hha', 0, 0, 'train', 'val'); 
  res_tr{end+1} = rcnn_all('task-detection', 'rgb_hha', 1, 0, 'train', 'val'); 
  res_tr{end+1} = rcnn_all('task-detection', 'rgb_hha', 1, 1, 'train', 'val'); 
  
  res{1} = rcnn_all('task-detection', 'rgb_hha', 0, 0, 'trainval', 'test');
  res{end+1} = rcnn_all('task-detection', 'rgb_hha', 1, 0, 'trainval', 'test');
  res{end+1} = rcnn_all('task-detection', 'rgb_hha', 1, 1, 'trainval', 'test');
  res{end+1} = rcnn_all('task-detection', 'rgb', 0, 0, 'trainval', 'test');
  res{end+1} = rcnn_all('task-detection', 'hha', 0, 0, 'trainval', 'test');
  global RUNNAME; RUNNAME = 'release';
  svmC = 10.^[-3 -4 -5 -6];
  for i = 1:length(svmC),
    res_svm{i} = rcnn_all('task-detection', 'rgb_hha', 0, 0, 'train', 'val', ...
      sprintf('_svmC%d', log10(svmC(i))), 'svm_C', svmC(i)); 
  end
  global RUNNAME; RUNNAME = 'release_ft-trainval';
  rcnn_all('task-detection', 'rgb_hha', 0, 0, 'test1', 'test2', '');
  rcnn_all('task-detection', 'rgb_hha', 0, 0, 'trainval', {'test2', 'test'}, 'matlab', '_matlab');
  rcnn_all('task-detection', 'rgb_hha', 1, 0, 'trainval', {'test2', 'test'}, 'matlab', '_matlab');
  rcnn_all('task-detection', 'rgb_hha', 1, 1, 'trainval', {'test2', 'test'}, 'matlab', '_matlab');
  rcnn_all('task-detection', 'rgb_hha', 0, 0, 'train', 'val', 'matlab', '_matlab');
  rcnn_all('task-detection', 'hha', 0, 0, 'test1', 'test2', 'matlab', '_matlab');
  rcnn_all('task-detection', 'rgb', 0, 0, 'test1', 'test2', 'matlab', '_matlab');
  rcnn_all('task-detection', 'rgb_hha', 0, 0, 'trainval', 'test', '_svmC-4', 'svm_C', 1e-4);
  
  global RUNNAME; RUNNAME = 'release';
  rcnn_all('task-detection', 'rgb_hha', 0, 0, 'test1', 'test2', 'release', '');
  rcnn_all('task-detection', 'hha', 0, 0, 'test1', 'test2', 'release', '');
  rcnn_all('task-detection', 'rgb', 0, 0, 'test1', 'test2', 'release', '');
end

if strcmp(jobName, 'region_vis')
  p = get_paths(); c = benchmarkPaths();
  NYU_ROOT_DIR = c.dataDir;
  REGIONDIR = fullfile(p.output_dir, 'regions', 'release-gt-inst');
  SALT = 'release'; MAX_BOXES = 2000;
  task = 'task-detection-with-cabinet';
  imset = 'val';
  
  imdb = imdb_from_nyud2(NYU_ROOT_DIR, imset, task, REGIONDIR, SALT, MAX_BOXES);
  imdb.roidb_func = @roidb_from_nyud2_region;
  roidb = imdb.roidb_func(imdb);
  box_dir = 'cache/release/detection/detector/rgb_hha_30000_train_region-features-1_region-task-1/detections/'; 
  out_dir = 'cache/release/detection/detector/rgb_hha_30000_train_region-features-1_region-task-1/vis/';
  cls = 'chair';
  for i = 1:19,
    vis_region_detect(imdb.classes{i}, box_dir, imdb, roidb, 0, 1, out_dir);
  end

end


if strcmp(jobName, 'compute_overlap')
  imlist = getImageSet('all');
  regiondir = 'cache/release/output/regions/release-gt-inst/';
  for i = 1:length(imlist),
    dt = load(fullfile_ext(regiondir, imlist{i}, 'mat'), 'bboxes', 'superpixels', 'sp2reg');
    [iu, inter, reg_area_1, reg_area_2] = compute_region_overlap(dt.superpixels, dt.sp2reg, dt.sp2reg);
    iu = sparse(double(iu)); save(fullfile_ext(regiondir, [imlist{i} '_iu'], 'mat'), '-append', 'iu');
  end
end
