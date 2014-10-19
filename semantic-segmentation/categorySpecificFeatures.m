function features = categorySpecificFeatures(imName, paths, param, models)
	load(fullfile(paths.ss_feature_dir, 'sift', strcat(imName, '.mat')), 'clusters');
	% [clusters] = getAmodalCompletion(imName);
	typ = param.featureParam.featureCacheName;

	for i = 1:size(clusters, 2),
		param.featureParam.selThresh = {i, i};
		[feature, superpixels, spArea] = loadExamples({imName}, paths, param, false);
		feature = feature{1};
		superpixels = superpixels{1};
		
		score = zeros(length(models), size(feature,2));
		for j = 1:length(models)
			score(j,:) = svmMulticlassTest(models(j), feature);
		end
		f = score;
		
		[cId, ind] = unique(clusters(:,i));
		features{i}(:, cId) = f(:, ind);
	end
	
	fileName = fullfile(paths.ss_feature_dir, typ, strcat(imName, '.mat'));
	save(fileName, 'clusters', 'superpixels', 'features');
end
