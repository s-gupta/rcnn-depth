function wrapperColorCues(inName, outName)
  I = imread(inName);

  dt = load('/work3/sgupta/eccv14-code/piotr-structured-edges-pami/models/forest/modelBsds.mat');
  model = dt.model;

  %% set detection parameters (can set after training)
  model.opts.multiscale=1;          % for top accuracy set multiscale=1
  model.opts.sharpen=0;             % for top speed set sharpen=0
  model.opts.nTreesEval=4;          % for top speed set nTreesEval=1
  model.opts.nThreads=4;            % max number threads for evaluation
  model.opts.nms=0;                 % set to true to enable nms


  [Es] = edgesDetect(I, model, zeros([size(I(:,:,1)), 0]));
  E = Es;
  save(outName, 'E', 'Es');
end
