function [E, Es, O] = detectEdge(Iin, Din, id, C, model, sc, outFile, cacheFile)
  if(isstr(Iin)), Iin = imread(Iin); end
  if(isstr(Din)), Din = imread(Din); end
  if(isstr(model)), load(model); end


  opts = model.opts;
  rgbd = opts.rgbd;

  model.opts.multiscale=1;          % for top accuracy set multiscale=1
  model.opts.nTreesEval=4;          % for top speed set nTreesEval=1
  model.opts.nThreads=1;            % max number threads for evaluation
  model.opts.sharpen=0;             % for top speed set sharpen=0

  if(rgbd == 3), 
    colorModel = model.colorModel.model;
    cues = getAllCues(Iin, Din, C, colorModel, opts.rgbd3opts.vars, cacheFile); 
  end
  
  for i = 1:length(sc),
    I = imresize(Iin, sc(i), 'lanczos3');
    ng = zeros(size(I,1), size(I,2), 0);

    if(rgbd), 
      D = imresize(Din, sc(i), 'nearest');
      D=single(D)/1e4;
    end
    if(rgbd==1), 
      I=D; 
    elseif(rgbd==2), 
      I=cat(3,single(I)/255,D);
    elseif(rgbd==3), 
      I=cat(3,single(I)/255,D,1e-3./D); 
      ng = cat(3, cues{:});
      ng = imresize(ng, sc(i), 'lanczos3');
    end
    model.opts.nms = 0;
    [Es{i}, O{i}] = edgesDetect(I,model,ng); 
    E{i} = edgesNmsMex(Es{i}, O{i}, 1, 5, 1.01, model.opts.nThreads);
    
    % model.opts.nms = 1;
    % [E{i} O{i}] = edgesDetect(I,model,ng); 
  end
  if(~isempty(outFile)), save(outFile, 'E', 'Es', 'O'); end
end
