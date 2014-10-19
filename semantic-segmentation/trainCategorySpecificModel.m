function trainCategorySpecificModel(imSet, paths, param, II)
	featureParam = param.featureParam;
	gtParam = param.gtParam;
	classifierParam = param.classifierParam;

	for i = 1:2,
		imList{i} = getImageSet(imSet{i});
	end

	%Load the features.
	cacheFile = fullfile(paths.ss_feature_dir, strcat(featureParam.featureCacheName, '_', param.fileSuffix, '.mat')); 
	if(exist(cacheFile, 'file') & param.useCache)
		dt = load(cacheFile);
		assert(isequal(imList, dt.imList), 'Image sets not same.... re-gathering features.. \n');
		feature = dt.feature; gt = dt.gt; sp = dt.sp;
		spArea = dt.spArea;
	else
		for i = 1:length(imList),
			[feature{i} sp{i} spArea{i} gt{i}] = loadExamples(imList{i}, paths, param, true);
		end
		save(cacheFile,'-v7.3','gt','feature','sp','spArea','imList', 'param', 'imSet');
	end
	fprintf('Loaded Features!\n');

	for i = 1:length(imList),
		gt{i} = cat(2, gt{i}{:});
		spArea{i} = cat(2, spArea{i}{:});
		feature{i} = cat(2, feature{i}{:});
	end

	fprintf('Preprocessing done!\n');
	

	%Train the SVMs here!
	for I = II,
		for i = 1:length(imList),
			gtLabel{i} = zeros(size(gt{i}));
			gtLabel{i}(gt{i} == I) = 1;
			gtLabel{i}(gt{i} ~= I) = 2;
			w{i} = spArea{i}; %ones(size(gt{i}));
		end

		%Train the model on train set
		modelFile = fullfile(paths.ss_model_dir, featureParam.featureCacheName, strcat(sprintf('%s_%02d', param.classifierFileName, I), '.mat')); 
		try
			load(modelFile)
		catch
			model = svmMulticlassTrain(feature, gtLabel, w, classifierParam);
			ap = model.ap(model.svmParamInd);
			save(modelFile, 'imList','param','model','ap','imSet');
		end
	end
end
