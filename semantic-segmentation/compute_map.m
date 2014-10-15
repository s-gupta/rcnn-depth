function [map, dimensions] = compute_map(imname, out_dir, param, data)
  try
    assert(~isempty(imname) && ~isempty(out_dir));
    fileName = fullfile_ext(out_dir, imname, 'mat');
    load(fileName, 'map', 'dimensions');
  catch
    switch(param.typ)
      case 'colorSift',
        % Call the binary to compute color sift here, for the time being loading features from disk here
        % dt = load(fullfile('/work4/sgupta/kinect/iResultsNYU/bowmaps/siftK1000_2000/', strcat(imname, '.mat')));
        % map = dt.bow_map;
        [map] = sift(data.I, param.siftParam, data.vocab); 

      case 'gTexton',
        % f we want to compute the geocentric textons then compute them here
        [map] = gTextons(data.pc, data.N, data.yDir, param);

    end
    dimensions = param.dimensions;

    if(~isempty(imname)),
      fileName = fullfile(out_dir, strcat(imname, '.mat'));
      save(fileName, 'map', 'dimensions');
    end
  end
end
