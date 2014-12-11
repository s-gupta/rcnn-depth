function [ov_cls, num_boxes] = benchmark_boxes(imdb, roidb, K, nms_thresh, plotstr) 
  for j = 1:length(imdb.classes), ov_cls{j} = []; end
  num_boxes = zeros(length(roidb.rois), length(K));
  for i = 1:length(roidb.rois),
    r = roidb.rois(i);
    boxes = r.boxes(r.gt == 0, :);
    
    if(nms_thresh < 1), 
      top = nmsOverlap(bboxOverlap(boxes, boxes), [size(boxes,1):-1:1]', nms_thresh);
      boxes = boxes(top, :);
    end
    ov = bboxOverlap(r.boxes(r.gt == 1,:), boxes);
    max_ov = zeros(size(ov,1), length(K)); 
    for k = 1:length(K),
      num_boxes(i,k) = min(K(k), size(ov,2));
      max_ov(:,k) = max(ov(:,1:num_boxes(i,k)), [], 2); 
    end
    cls = r.class(r.gt == 1); 
    
    for j = 1:length(cls),
      if(cls(j) > 0), ov_cls{cls(j)}(end+1,:) = max_ov(j,:); end 
    end
    fprintf('.');
  end
  num_boxes = mean(num_boxes, 1);

  for j = 1:length(imdb.classes),
    subplot(4,5,j); 
    plot(num_boxes, mean(ov_cls{j} > 0.7), plotstr);
    grid on; grid minor;
    % xlabel('Number of candidates', 'FontSize', 18);
    % ylabel('Fraction of gt instance with more than x overlap', 'FontSize', 18);
    axis([min(K) max(K), 0, 1.0]);
    set(gca,'XScale','log');
    % set(gca, 'FontSize', 18);
    title(sprintf('%s', imdb.classes{j})); 
  end
end
