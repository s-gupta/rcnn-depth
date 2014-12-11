function startup()
  % genPath = {'piotr-structured-edges-pami', 'rgbdutils', 'utils'};
  if(~exist('nyu-hooks', 'dir')), warning('Missing nyu-hooks. You will need it for some functions. Get it from here: https://github.com/s-gupta/nyu-hooks'); end
  if(~exist('utils', 'dir')), warning('Missing utils. You will need it for some fuctions. Get it from here: https://github.com/s-gupta/utils/'); end
  if(~exist('rgbdutils', 'dir')), warning('Missing rgbd-utils. You will need it for some functions. Get it from here: https://github.com/s-gupta/rgbdutils'); end

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
