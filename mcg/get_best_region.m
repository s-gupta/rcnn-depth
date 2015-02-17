function [IIgt, II] = get_best_region(imName, oid)
%%
% outDir = 'rmrc-talk';
% imName = 'img_6449';
% for i = 1:20,
%  [IIgt, II] = get_best_region(imName, i);
%  imwrite(IIgt, fullfile_ext(outDir, sprintf('%s_%02d_gt', imName, i), 'png'));
%  imwrite(II, fullfile_ext(outDir, sprintf('%s_%02d_pr', imName, i), 'png'));
% end
  REGIONDIR = '~/psi/work4/sgupta/regions-nyu/datasets/nyud40Obj/mcg_sfdnV03_3_3_depthV2_030/';
  gt = getGroundTruth(imName, 'instance');
  I = getImage(imName, 'images');
  mask = gt == oid;
  dt = load(fullfile(REGIONDIR, [imName '_cands.mat']));
  for i = 1:500, %size(sp2reg,1),
    spi = dt.sp2reg(:,i);
    maskReg = spi(dt.sp);
    ii(i) = nnz(maskReg(:) & mask(:));
    uu(i) = nnz(maskReg(:) | mask(:));
  end
  o = ii./uu;
  [~, ind] = max(o);
  spi = dt.sp2reg(:,ind);
  maskReg = spi(dt.sp);
  figure(1); 
  subplot(2,2,1); imagesc(maskReg);
  IIgt = im2double(overlay(im2double(I), mask, zeros(size(mask))));
  [ii, jj ] = find(mask);
  bbox(1) = min(jj); bbox(2) = min(ii); bbox(3) = max(jj); bbox(4) = max(ii);
  IIgt = draw_rect_vec(IIgt, bbox', [0 1 0], 3);

  subplot(2,2,2); imagesc(IIgt);
  % II = im2double(overlay(repmat(rgb2gray(im2double(I)), [1 1 3]), maskReg, zeros(size(mask))));
  II = im2double(overlay(im2double(I), maskReg, zeros(size(mask))));
  II = draw_rect_vec(II, dt.bboxes(ind, [2 1 4 3])', [0 1 0], 3);
  subplot(2,2,3); imagesc(II);
  subplot(2,2,4); imagesc(gt);
end

