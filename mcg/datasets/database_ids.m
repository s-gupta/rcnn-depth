function ids = database_ids( database, gt_set )
    switch database,
      case 'nyud4oObj',
        ids = getImageSet(gt_set);
      otherwise,
        index_file = fullfile(mcg_root_dir,'datasets',database, 'gt_sets',[gt_set '.txt']);
        fileID = fopen(index_file);
        ids = textscan(fileID, '%s');
        ids = ids{1};
        fclose(fileID);
    end
end

