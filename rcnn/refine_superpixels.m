function candidates_new = refine_superpixels(inst, candidates, N)
% function candidates_new = refine_superpixels(inst, candidates, N)

  sp=candidates.superpixels;
  N1=min(N, numel(candidates.labels));
  reg2sp=false(max(sp(:)),N1); % nSP x nR
  for j = 1:N1
    reg2sp(candidates.labels{j},j)=true;
  end
  old.superpixels = sp;
  old.sp2reg = reg2sp';
  old.labels = candidates.labels(1:N1);
  
  % compress to remove irrelevant sps
  [sp, reg2sp] = compressSp2reg(sp, reg2sp);


  % add ground truth to reg2sp
  num_inst = max(inst(:));
  gt_sp = double(inst+1);
  gt_reg2sp = [zeros(1, num_inst); eye(num_inst)];
  gt_reg2sp = logical(gt_reg2sp);
  [sp, reg2sp] = intersect_partitions(gt_sp, sp, gt_reg2sp, reg2sp);

  candidates_new.scores = candidates.scores(1:N1,:);
  candidates_new.bboxes = candidates.bboxes(1:N1,:);
  candidates_new.superpixels = sp;
  candidates_new.sp2reg = reg2sp';
  candidates_new.sp2reg = candidates_new.sp2reg(num_inst+1:end,:);
  candidates_new.labels = cell(N1,1);
  candidates_new.old = old; 
  for i = 1:N1,
    candidates_new.labels{i} = uint32(find(candidates_new.sp2reg(i,:)));
  end

  %% Test that the regions are the same still
  for i = 1:N1,
    sp2reg1 = false(1, max(candidates.superpixels(:)));
    sp2reg1(candidates.labels{i}) = true;

    sp2reg2 = candidates_new.sp2reg(i,:);
    % figure(1); subplot(1,2,1); imagesc(sp2reg2(candidates_new.superpixels));
    % figure(1); subplot(1,2,2); imagesc(sp2reg1(candidates.superpixels));
    assert(isequal(sp2reg2(candidates_new.superpixels), sp2reg1(candidates.superpixels)), sprintf('%d not good.', i));
  end
end
      
function [sp, sp2reg] = intersect_partitions(sp1, sp2, sp2reg1, sp2reg2)
  %do a unique
  [unique_pairs, junk, sp]=unique([sp1(:) sp2(:)], 'rows');
  sp=reshape(sp, size(sp1));

  %create new sp2reg
  sp2reg=[sp2reg1(unique_pairs(:,1),:) sp2reg2(unique_pairs(:,2),:)];

  %break disconnected sp
  spNew = zeros(size(sp));
  sp2regNew = false(size(sp2reg));
  cnt = 0;
  %% Check connectivity of superpixels
  for i = 1:max(sp(:))
    tt = bwconncomp(sp == i);
    for j = 1:tt.NumObjects,
      cnt = cnt+1;
      spNew(tt.PixelIdxList{j}) = cnt;
      sp2regNew(cnt, :) = sp2reg(i, :);
    end
  end

  sp=spNew;
  sp2reg=sp2regNew;
end

function [spNew sp2regNew] = compressSp2reg(sp, sp2reg)
  assert(islogical(sp2reg), 'sp2reg should be logical');
  [sp2reg bb spSame] = unique(sp2reg, 'rows');
  sp = spSame(sp);

  spNew = zeros(size(sp));
  sp2regNew = false(size(sp2reg));
  cnt = 0;
  %% Check connectivity of superpixels
  for i = 1:max(sp(:))
    tt = bwconncomp(sp == i);
    for j = 1:tt.NumObjects,
      cnt = cnt+1;
      spNew(tt.PixelIdxList{j}) = cnt;
      sp2regNew(cnt, :) = sp2reg(i, :);
    end
  end
end 
