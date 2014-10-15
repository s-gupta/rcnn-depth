paths = getPaths();

%% Create the file for storing the results
allResultsFileName = fullfile(paths.outDir, 'allSSResults.mat');
if(~exist(allResultsFileName, 'file'))
	save(allResultsFileName, 'allResultsFileName');
end

%% Get the ucmThresh
%  dt = load(fullfile(paths.modelDir, 'allBUSResults.mat'));
%  ucmThresh = dt.th.ucmThresh;

f = struct('computeGeneric', true, 'computeCategorySpecific', true, 'ablationStudy', true, ...
	'full', true, 'scene', true, 'sceneF', true, 'fullScene', true);

% Compute the genericlassifier, sift and g texton features
if(f.computeGeneric)
	wrapperComputeFeatures1('train', ucmThresh);
	wrapperComputeFeatures1('val', ucmThresh);
	wrapperComputeFeatures1('test', ucmThresh);
end

% Compute the category specific feaatures - requries some training.
if(f.computeCategorySpecific)
	typ = {'categorySpecific-level1-all', 'categorySpecific-level1-gTexton' ,'categorySpecific-level1-sift'},
	for i = 1:length(typ), 
		wrapperComputeFeatures2(typ{i}, 'train');
		wrapperComputeFeatures2(typ{i}, 'val');
		wrapperComputeFeatures2(typ{i}, 'test');
	end
end


% Train the models and generate numbers in the paper
if(f.ablationStudy)
	classMapping = 'classMapping40';
	trSet = 'train1';
	valSet = 'train2';
	testSet = 'val';
	classifier = {'svm-full', 'tree-full'},
	typ = {'full', 'generic', 'categorySpecific', 'onlyGeom', 'onlyApp', 'noAmodal'},

	for i = 1:2,
		for j = 1:6, 
			evalResAblation(i,j) = wrapperTrainTestBenchmarkModel(trSet, valSet, testSet, typ{j}, classifier{i}, classMapping);
			fprintf('%s - %s, accuracy found to be, %0.3f.\n', classifier{i}, typ{j}, evalResAblation(i,j).fwavacc);
		end
	end
	save(allResultsFileName, '-APPEND', 'evalResAblation');
end


%% Train the full models on the whole data here
if(f.full)
	classMapping = 'classMapping40';
	trSet = 'train';
	valSet = 'val';
	testSet = 'test';
	typ = 'full';
	classMapping = {'classMapping40', 'classMapping04'};
	classifier = {'svm-full', 'tree-full'};
	for i = 1:2,
		for j = 1:2,
			evalResFull(i,j) = wrapperTrainTestBenchmarkModel(trSet, valSet, testSet, typ, classifier{i}, classMapping{j});
			fprintf('%s - %s, accuracy found to be, %0.3f.\n', classifier{i}, classMapping{j}, evalResFull(i,j).fwavacc);
		end
	end
	save(allResultsFileName, '-APPEND', 'evalResFull');
end

%% Run the secne classification stuff
if(f.scene)
	modelFileName = fullfile(paths.modelDir, 'svm-full_entryLevel_ancFull_tr-train_val-val_useVal-1');
	[softOutputDir, hardOutputDir] = testModel('train', paths, modelFileName);
	[softOutputDir, hardOutputDir] = testModel('val', paths, modelFileName);
	[softOutputDir, hardOutputDir] = testModel('test', paths, modelFileName);
	objectDir = softOutputDir;
	trSet = 'train';
	valSet = 'val';
	testSet = 'test';
	sceneMapping = 'sceneMapping10';
	classMapping = 'classMapping40';
	typ = {'gTexton', 'colorSift', 'gTextonSift', 'objScores'};

	for i = 1:4,
		evalResScene(i) = wrapperScene(typ{i}, trSet, valSet, testSet, classMapping, sceneMapping, objectDir);
		fprintf('%s - scene caltechAccuracy: %0.3f.\n', typ{i}, evalResScene(i).caltechAccuracy);
	end
	save(allResultsFileName, '-APPEND', 'evalResScene');
end

%% Run the scene classifiers to compute the scene scores 
if(f.sceneF)
	sceneModelName = fullfile(paths.modelDir, 'scene-objScores_tr-train_val-val_useVal-1'); %TODO
	outputFileName = testSceneModel('train', paths, sceneModelName);
	wrapperComputeFeatures3('train', outputFileName);
	outputFileName = testSceneModel('val', paths, sceneModelName);
	wrapperComputeFeatures3('val', outputFileName);
	outputFileName = testSceneModel('test', paths, sceneModelName);
	wrapperComputeFeatures3('test', outputFileName);
end

%% Train the full + scene models 
if(f.fullScene)
	classMapping = 'classMapping40';
	trSet = 'train';
	valSet = 'val';
	testSet = 'test';
	typ = 'full+scene';
	classifier = {'svm-full', 'tree-full'},
	classMapping = {'classMapping40', 'classMapping04'},
	for i = 1:2,
		for j = 1:2,
			evalResFullScene(i,j) = wrapperTrainTestBenchmarkModel(trSet, valSet, testSet, typ, classifier{i}, classMapping{j});
			fprintf('%s - %s, accuracy found to be, %0.3f.\n', classifier{i}, classMapping{j}, evalResFullScene(i,j).fwavacc);
		end
	end
	save(allResultsFileName, '-APPEND', 'evalResFullScene');
end
