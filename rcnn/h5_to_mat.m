function h5_to_mat(h5_file, imdb, window_file, output_dir)
  pt = load(window_file);
  roidb = imdb.roidb_func(imdb);
  
  a_counter = 0;
  a = read_h5_file(h5_file, sprintf('data-%06d', a_counter));
  a = a{1};
  rem = size(a, 4);
  for i = 1:length(roidb.rois),
    d = roidb.rois(i);
    assert(isequal(pt.imlist{i}, imdb.image_ids{i}));
    assert(isequal(single(pt.list{i}(:,3:6)+1), d.boxes));
   
    while rem < size(d.boxes,1),
      a_counter = a_counter+1;
      a_new = read_h5_file(h5_file, sprintf('data-%06d', a_counter));
      a = cat(4, a, a_new{1});
      rem = size(a,4);
      fprintf('%d (%d), ', rem, a_counter);
    end
    fprintf('\n');
    f = a(:,:,:,[1:size(d.boxes,1)]);
    a = a(:,:,:, size(d.boxes,1)+1:end);
    rem = size(a,4);
    d.feat = permute(f, [4 3 1 2]);
    save_file = [output_dir imdb.image_ids{i} '.mat'];
    save(save_file, '-STRUCT', 'd');
  end
end
