function wrapper_rgbd_to_ucm(imset, modelfile, ucm_dir, contours_cues_dir, C)
  sc = [2 1 0.5];
  imlist = getImageSet(imset);
  
  model = load(modelfile); %fullfile_ext(contours_model_dir, 'forest', modelname, 'mat')); 
  model = model.model;
  
  exists_or_mkdir(fullfile(ucm_dir, 'multi'));
  exists_or_mkdir(fullfile(ucm_dir, 'scale_2.0'));
  exists_or_mkdir(fullfile(ucm_dir, 'scale_1.0'));
  exists_or_mkdir(fullfile(ucm_dir, 'edges_1.0'));
  exists_or_mkdir(fullfile(ucm_dir, 'scale_0.5'));
  exists_or_mkdir(fullfile(ucm_dir, 'multi-png'));

  for i = 1:length(imlist),
    id = imlist{i};
    I = getImage(id, 'images');
    D = getImage(id, 'depth');
    cacheFile = fullfile_ext(contours_cues_dir, id, 'mat');
    [E, Es, O] = detectEdge(I, D, [], C, model, sc, [], cacheFile);
    [ucm2 ucms] = contours_to_ucm(I, sc, Es, O);
    parsave(fullfile_ext(ucm_dir, 'multi', id, 'mat'), 'ucm2', ucm2);
    parsave(fullfile_ext(ucm_dir, 'edges_1.0', id, 'mat'), 'E', E{2}, 'Es', Es{2}, 'O', O{2});
    parsave(fullfile_ext(ucm_dir, 'scale_2.0', id, 'mat'), 'ucm2', ucms(:,:,1));
    parsave(fullfile_ext(ucm_dir, 'scale_1.0', id, 'mat'), 'ucm2', ucms(:,:,2));
    parsave(fullfile_ext(ucm_dir, 'scale_0.5', id, 'mat'), 'ucm2', ucms(:,:,3));
    imwrite(im2uint8(ucm2(3:2:end, 3:2:end)), fullfile_ext(ucm_dir, 'multi-png', id, 'png'));
  end 
end
