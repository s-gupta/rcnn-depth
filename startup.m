function startup()
  % genPath = {'piotr-structured-edges-pami', 'rgbdutils', 'utils'};
  addPath = {'nyu-hooks', 'scripts', 'mcg', 'utils', 'caffe/matlab/caffe/', 'rgbdutils'};
  genPath = {'structured-edges', 'rcnn'};

  for i = 1:length(addPath),
    addpath(addPath{i});
  end

  for i = 1:length(genPath),
    addpath(genpath_exclude(genPath{i}, '.git'));
  end
  startup_utils;
  startup_rgbdutils;
  startup_mcg;
  fprintf('startup done\n');
end
