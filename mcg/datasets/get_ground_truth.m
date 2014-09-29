function ground_truth = get_ground_truth( database, image_id )
    if strcmp(database,'pascal2012')
        ground_truth.object = imread(fullfile(database_root_dir(database), 'SegmentationObject', [image_id '.png']));
        ground_truth.class  = imread(fullfile(database_root_dir(database), 'SegmentationClass', [image_id '.png']));
    elseif strcmp(database,'bsds500')
        ground_truth = loadvar(fullfile(database_root_dir(database), 'ground_truth', [image_id '.mat']),'gt_seg');
    elseif strcmp(database,'nyud40Obj')
        ground_truth.object = getGroundTruth(image_id, 'instance');
        ground_truth.class  = getGroundTruth(image_id, 'segmentation', 'classMapping40');
        ground_truth.class(ismember(ground_truth.class, [1 2 8 22 38])) = 0; %tmp.groundTruth{1}.SegmentationClass;
        % Remove class 0 instances
        ground_truth.object(ground_truth.class == 0) = 0;
        obj = ground_truth.object;
        newobj = zeros(size(obj));
        uq = unique(obj(obj > 0));
        for i = 1:length(uq),
          newobj(obj == uq(i)) = i;
        end
        ground_truth.object = newobj;
    else
        error(['Unknown database: ' database]);
    end
end

