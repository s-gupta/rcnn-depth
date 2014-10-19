function [softOutputDir hardOutputDir] = testModel(imSet, paths, modelFileName)
	dt = load(modelFileName);
	param = dt.param;
	
	featureParam = param.featureParam;
	gtParam = param.gtParam;
	classifierParam = param.classifierParam;
	
	imList = getImageSet(imSet);

	[feature, sp, spArea] = loadExamples(imList, paths, struct('featureParam', featureParam, 'gtParam', gtParam), false);

	softOutputDir = fullfile(paths.ss_dir, 'soft_outs', strcat(sprintf('%s_%s', param.classifierFileName, param.fileSuffix)));
	hardOutputDir = fullfile(paths.ss_dir, 'hard_outs', strcat(sprintf('%s_%s', param.classifierFileName, param.fileSuffix)));
	
	if(~exist(softOutputDir, 'dir')) mkdir(softOutputDir); end
	if(~exist(hardOutputDir, 'dir')) mkdir(hardOutputDir); end

	parfor i = 1:length(imList),
		f = feature{i};
		s2 = []; p2 = [];
	
		switch param.classifierType,
			case 'svm-full',
				[s2, p2, ~] = svmMulticlassTest(dt.model, f);	
			case 'tree-full',
				p2 = testTreeBagging(dt.model, f);
				s2 = p2;
		end

		[score cls] = max(p2,[],1);
		superpixels = sp{i};
		scores = p2;
		rawScores = s2; 

		fileName = fullfile(softOutputDir, strcat(imList{i}, '.mat'));
		parsave(fileName, 'superpixels', superpixels, 'scores', scores, 'rawScores', rawScores);

		if(isfield(dt, 'thresh'))
			% Write to the output directory!
			spLabel = cls .* (score >= dt.thresh);
			segmentation = spLabel(sp{i});
			fileName = fullfile(hardOutputDir, strcat(imList{i}, '.mat'));
			parsave(fileName, 'segmentation', segmentation);
		end
	end
end
