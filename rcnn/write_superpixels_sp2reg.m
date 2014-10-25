function write_superpixels_sp2reg(imdb, roidb, spdir, sp2regdir)
% function write_superpixels_sp2reg(imdb, roidb, spdir, sp2regdir)

  rois = roidb.rois; 
  for i = 1:length(rois),
    assert(max(rois(i).sp(:)) < 10000);  
    imwrite(uint16(rois(i).sp), fullfile_ext(spdir, imdb.image_ids{i}, 'png'));
    imwrite(im2uint8(rois(i).sp2reg), fullfile_ext(sp2regdir, imdb.image_ids{i}, 'png'));
  end
end
