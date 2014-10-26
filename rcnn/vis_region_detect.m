function vis_region_detect(cls, box_dir, imdb, roidb, vis, saveit, out_dir)
  nms_thresh = 0;
  [~, clsId] = ismember(cls, imdb.classes);

  mkdir_if_missing(fullfile(out_dir, imdb.imset, cls));
  out_dir = fullfile(out_dir, imdb.imset, cls);
 
  dt = load(fullfile(box_dir, sprintf('%s_boxes_%s.mat', cls, imdb.name)));
  boxes = dt.boxes;
  rois = roidb.rois;
  % Load the ground truth structures, for the imdb
  parfor i = 1:length(imdb.image_ids),
    roi = rois(i);

    % Do non max suppression
    sp2regI = roi.sp2reg(roi.gt == 0, :);
    bboxI = roi.boxes(roi.gt == 0, :);
    [iu, ~, ~, ~] = compute_region_overlap(roi.sp, sp2regI, sp2regI);
    assert(isequal(boxes{i}(:,1:4), bboxI));
    scI = boxes{i}(:,end);
    pick = nmsOverlap(iu, scI, nms_thresh);
    
    sc{i} = cat(2, scI(pick), pick);
    sc{i}(:,3) = i;
  end

  sc = cat(1, sc{:});
  [~, ind] = sortrows(sc, -1);
  sc = sc(ind,:);

  for i = 1:200, %size(sc,1),
    image_id = sc(i,3); box_id = sc(i,2);
    roi = roidb.rois(image_id);
    I = getImage(imdb.image_ids{image_id}, 'images');
    sp2reg = roi.sp2reg(roi.gt == 0, :);
    sp = roi.sp;
    boxes = roi.boxes(roi.gt == 0, :);

    
    sp2regi = sp2reg(box_id,:);
    mask = sp2regi(sp);
    II = im2uint8(drawRegionsPaper(I, mask, 1));
    II = draw_rect_vec(II, boxes(box_id,:)', [0 255 255], 1);
    
    if vis,
      figure(1); 
      imagesc(cat(2, I, II)); axis image;
      % subplot(1,2,1); imagesc(I); plotDS(boxes(box_id,:), 'r'); axis image;
      % subplot(1,2,2); sp2regi = sp2reg(box_id, :); imagesc(mask); axis image
      pause;
    end
    if saveit,
      imwrite(cat(2, I, II), fullfile(out_dir, sprintf('vis_%03d.jpg', i)));
    end
  end
end
