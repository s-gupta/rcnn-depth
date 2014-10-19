function modelFile = trainModel(imSet, paths, param)
  
  featureParam = param.featureParam;
  gtParam = param.gtParam;
  classifierParam = param.classifierParam;
  bParam = param.bParam;

  for i = 1:2,
    imList{i} = getImageSet(imSet{i});
  end

  % Loading the features here
  cacheFile = fullfile(paths.ss_feature_dir, strcat(featureParam.featureCacheName, '_', param.fileSuffix, '.mat'))
  if(exist(cacheFile, 'file') & param.useCache)
    dt = load(cacheFile);
    assert(isequal(imList, dt.imList), 'Image sets not same.... re-gathering features.. \n');
    feature = dt.feature; gt = dt.gt; sp = dt.sp;
    spArea = dt.spArea;
  else
    for i = 1:length(imList),
      [feature{i} sp{i} spArea{i} gt{i}] = loadExamples(imList{i}, paths, struct('featureParam', featureParam, 'gtParam', gtParam), true);
    end
    save(cacheFile,'-v7.3','gt','feature','sp','spArea','imList','param','imSet');
  end
  fprintf('Loaded Features!\n');

  for i = 1:length(imList),
    gt{i} = cat(2, gt{i}{:});
    spArea{i} = cat(2, spArea{i}{:});
    feature{i} = cat(2, feature{i}{:});
  end

  % Training the model here
  modelFile = fullfile(paths.ss_model_dir, strcat(sprintf('%s_%s', param.classifierFileName, param.fileSuffix), '.mat')); 
  try
    dt = load(modelFile);
    model = dt.model;
  catch
    switch param.classifierType
        case 'svm-full',
          model = svmMulticlassTrain(feature, gt, spArea, classifierParam);
        case 'tree-full',
          model = trainTreeBagging(feature, gt, spArea, classifierParam);
    end
     save(modelFile, 'model', 'param', 'imList', 'imSet');
  end

  %% Test the model on the validation data here and pick the threshold for projection and store it again in the model
  try
    dt = load(modelFile);
    thresh = dt.thresh;
  catch
    outputDir = testModel(imSet{2}, paths, modelFile);
    
    evalResVal = benchmarkSSInternal(outputDir, imSet{2}, bParam);
    [~, ind] = max(evalResVal.detailsIU.fwavacc);
    thresh = bParam.threshIU(ind);
    fprintf('Threshold found to be %0.2f\n', thresh);
    save(modelFile, '-append', 'thresh', 'evalResVal');
  end
end
