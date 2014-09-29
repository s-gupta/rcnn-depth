function cues = getAllCues(Iin, Din, C, colorModel, vars, cacheFile)
  cache = ~isempty(cacheFile);
  try
    assert(cache);
    load(cacheFile, vars{:});
  catch
    D = double(Din)./10;
    [y1, y2, y3, angl1, angl2] = yCues(D, C, 1);
    [ng1, ng2, dg] = normalCues(D, C, 1);

    model = colorModel;
    %% set detection parameters (can set after training)
    model.opts.multiscale=1;          % for top accuracy set multiscale=1
    model.opts.sharpen=0;             % for top speed set sharpen=0
    model.opts.nTreesEval=4;          % for top speed set nTreesEval=1
    model.opts.nThreads=4;            % max number threads for evaluation
    model.opts.nms=0;                 % set to true to enable nms

    Es = edgesDetect(Iin, model, zeros([size(Iin(:,:,1)), 0]));
    if(cache), save(cacheFile, 'y1', 'y2', 'y3', 'angl1', 'angl2', 'ng1', 'ng2', 'dg', 'Es'); end
  end
   
  for ii = 1:length(vars),
    eval(sprintf('cues{ii} = %s;', vars{ii}));
    %% Add the ng at scale 1 and 2 and use that for scale 3
    if(size(cues{ii},4) > 1), cues{ii}(:,:,:,3) = cues{ii}(:,:,:,1)+cues{ii}(:,:,:,2); end
    cues{ii}(:,:,end+1,:) = max(cues{ii}, [], 3);
    cues{ii} = reshape(cues{ii}, size(cues{ii},1), size(cues{ii},2), []);
  end
end
