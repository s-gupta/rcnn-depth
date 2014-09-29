function p = get_paths(runname)
  if ~exist('runname', 'var')
    global RUNNAME
    if isempty(RUNNAME), runname = 'release'
    else, runname = RUNNAME; end
  end
  p.root_cache_dir = fullfile(pwd(), '..', 'eccv14-cachedir');
  p.cache_dir = fullfile(p.root_cache_dir, runname);
  
  p.contours_dir = fullfile(p.cache_dir, 'contours');
    p.contours_cues_dir = fullfile(p.contours_dir, 'cues');
    p.contours_model_dir = fullfile(p.contours_dir, 'models');
  
  p.output_dir = fullfile(p.cache_dir, 'output');
    p.ucm_dir = fullfile(p.output_dir, 'ucm'); 
  p.box_dir = fullfile(p.cache_dir, 'regions');
  
  % p.caffe_dir = '/home/eecs/sgupta/eccv14-code/cachedir/release/detection/';
  %   p.ft_image_dir = fullfile(p.caffe_dir, 'ft_images');
  %   p.ft_hha_dir = fullfile(p.caffe_dir, 'ft_hha');
  
  p.detectior_dir = fullfile(p.cache_dir, 'detection');
  p.cnnF_cache_dir = fullfile(p.detectior_dir, 'feat_cache');
    p.ft_image_dir = fullfile(p.detectior_dir, 'ft_images');
    p.ft_hha_dir = fullfile(p.detectior_dir, 'ft_hha');
    p.ft_disparity_dir = fullfile(p.detectior_dir, 'ft_disparity');
    p.ft_dir = fullfile(p.detectior_dir, 'finetuning');
      p.wf_dir = fullfile(p.ft_dir, 'windowfile');
      p.proto_dir = fullfile(p.ft_dir, 'protodir');
      p.snapshot_dir = fullfile(p.ft_dir, 'snapshot');

  f = fieldnames(p);
  for i = 1:length(f), exists_or_mkdir(p.(f{i})); end
  
  p.mean_file_color = fullfile('caffe-data', 'mean', 'imagenet_mean');
  p.mean_file_disparity = fullfile('caffe-data', 'mean', 'nyu_disparity_train1');
  p.mean_file_hha = fullfile('caffe-data', 'mean', 'nyu_hha_train1');
  p.caffe_net = fullfile('caffe-data', 'caffe_reference_imagenet_model');
  p.rgb_edge_model = fullfile('structured-edges', 'models', 'forest', 'modelBsds.mat');
end
