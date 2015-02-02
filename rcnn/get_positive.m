% ------------------------------------------------------------------------
function [X_pos, keys] = get_positive(imdb, rcnn_model, do_latent)
% ------------------------------------------------------------------------
opts = rcnn_model.training_opts;

X_pos = cell(max(imdb.class_ids), 1);
keys = cell(max(imdb.class_ids), 1);

class_ids = imdb.class_ids;
for j = imdb.class_ids
  if isempty(X_pos{j})
    X_pos{j} = single([]);
    keys{j} = [];
  end
end
 
for i = 1:length(imdb.image_ids)
  tic_toc_print('%s (rcnn_region_train) : pos features %d/%d\n', ...
                procid(), i, length(imdb.image_ids));
  if do_latent 
    di = get_features(imdb, i, rcnn_model.feat_opts, false);
    
    % Code to select the examples according to overlap with ground truth...
    [iu, inter, reg_area_1, reg_area_2] = compute_region_overlap(di.sp, di.sp2reg, di.sp2reg(di.gt, :));
    iu = iu > opts.pos_ov_thresh;
    feat = rcnn_scale_features(di.feat, opts.feat_norm_mean);
    zs = bsxfun(@plus, feat*rcnn_model.detectors.W, rcnn_model.detectors.B);
    INFINF = 2*max(zs(:));
    % Code to pick the highest scoring detection which overlaps the largest with this ground truth instance.
    feats_sel = di.feat(di.gt, :);
    feats_sel(:) = NaN;
    for j = 1:size(iu,1),
      if(ismember(di.class(j), class_ids))
        cls = di.class(j);
        sc = (zs(:,di.class(j)) + INFINF) .* iu(:,j) - INFINF;
        [~, ind] = max(sc);
        feats_sel(j,:) = di.feat(ind,:);
        X_pos{cls} = cat(1, X_pos{cls}, di.feat(ind,:));
        keys{cls} = cat(1, keys{cls}, [i ind]);
      end
    end
  else
    di = get_features(imdb, i, rcnn_model.feat_opts, true);
    for j = imdb.class_ids
      sel = find(di.class == j);
      if ~isempty(sel)
        X_pos{j} = cat(1, X_pos{j}, di.feat(sel,:));
        keys{j} = cat(1, keys{j}, [i*ones(length(sel),1) sel]);
      end
    end
  end
end
