function [ex, ds, thresh, recall, ind] = load_detector(cls, detDir, threshSet, detSet, pThresh, nmsThresh, boxFieldName)
% function [ex, ds, thresh, recall, ind] = load_detector(cls, detDir, threshSet, detSet, pThresh, nmsThresh, boxFieldName)
%   nmsThresh - is the overlap used to nms suppress the output (for deep detector outputs)
%

% AUTORIGHTS

  try
    cls = regexprep(cls, ' ', '-');
    fprintf('load_detector: Loading pr curve from %s.\n', fullfile(detDir, sprintf('%s_boxes_nyud2_%s_release.mat', cls, threshSet)));
    dt1 = load(fullfile(detDir, sprintf('%s_boxes_nyud2_%s_release.mat', cls, threshSet)), boxFieldName);
    dt2 = load(fullfile(detDir, sprintf('%s_pr_nyud2_%s_release.mat', cls, threshSet)), 'prec', 'recall');
   
    dt2.prec = [1; dt2.prec];
    fprintf('load_detector: Doing nms for %s for %s.\n', cls, threshSet);
    ds = cellfun(@(x) x(rcnn_nms(x, nmsThresh), :), dt1.(boxFieldName), 'UniformOutput', false);
    ds = cat(1, ds{:});
    dt2.recall = [0; dt2.recall];
    ds = [0 0 0 0 inf; ds];

    %dt2.prec = [1; dt2.prec];
    %ds = cat(1, dt1.ds{:});

    ind = find(dt2.prec > pThresh, 1, 'last');
    recall = dt2.recall(ind);
    sc = sort(ds(:,5), 'descend');
    thresh = sc(ind);
    
    dt1 = load(fullfile(detDir, sprintf('%s_boxes_nyud2_%s_release.mat', cls, detSet)), boxFieldName);
    fprintf('load_detector: Doing nms for %s for %s.\n', cls, detSet);
    ds = cellfun(@(x) x(rcnn_nms(x, nmsThresh), :), dt1.(boxFieldName), 'UniformOutput', false);
    ex = true;
  
  catch e
    prettyexception(e);
    error('load_detector: Detector for %s does not exist!\n', cls);
    ds = []; thresh = inf; recall = 0; ex = false;
  end
end
