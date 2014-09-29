function image = get_image( database, image_id )
    if strcmp(database,'pascal2012')
        image = imread(fullfile(database_root_dir(database), 'JPEGImages', [image_id '.jpg']));
    elseif strcmp(database,'bsds500')
        image = imread(fullfile(database_root_dir(database), 'images', [image_id '.jpg']));
    else
        error(['Unknown database: ' database]);
    end
end

