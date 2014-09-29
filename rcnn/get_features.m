function dOut = get_features(imdb, image_ids_i, opts, onlyGt)
% function get_features(imName, opts)
% Various types of features, as encoded by opts

% AUTORIGHTS

  for i = 1:length(opts),
    d{i} = get_features_i(imdb, image_ids_i, opts(i), onlyGt);
    % if(isfield(opts(i).param, 'wt'))
    %   d{i}.feat = (opts(i).param.wt).*d{i}.feat;
    % end
  end
  for i = 2:length(d),
    assert(isequal(d{1}.boxes, d{i}.boxes));
  end
  dOut = d{1};
  for i = 1:length(d),
    dd(i).feat = d{i}.feat;
  end
  dOut.feat = cat(2, dd(:).feat);
  assert(any(isnan(dOut.feat(:))) == false, 'Features were NaN!!');
end
