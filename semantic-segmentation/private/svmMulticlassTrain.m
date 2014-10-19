function model = svmMulticlassTrain(X, Y, W, param, fileName, svmParam, mlrParam)
% function svmMulticlassTrain(X, y, param)
% Input: 
% param.iksvmN  - 
% param.C   - array of C to use
% param.w1  - array of w1 to use
% param.bias  - array of bias to use
% param.numClass - number of classes 1 to N
% param.CMLR
% param.biasMLR
% param.mapNanToZero
% param.useValForTraining
% X{1} is training and X{2} is validation set of size DxN
% Y{1} must be of size 1xN, W{1} must be of size 1x
  param

  if(~exist('fileName', 'var')) fileName = sprintf('%d-%s.mat', getPID(), datestr(now, 'yyyymmddTHHMMSS')); end

  for i = 1:2,
    assert(size(Y{i},1) == 1, 'Y must be a 1xN vector.');
    assert(size(Y{i},2) == size(X{i},2), 'Number of instances do not match number of labels supplied.');
    assert(size(Y{i},2) == size(W{i},2), 'Number of instances do not match number of weights supplied.');
  end 

  model.isSparse = param.iksvmN == 0;
  model.iksvmN = param.iksvmN;
  model.numClass = param.numClass;
  model.mapNanToZero = param.mapNanToZero;
  model.weighted = param.weighted;
  % Train SVMs on the training set.
  trVal = 1;
  [~, model.fNorm] = normalizeFeatures(X{trVal});
 
  if((exist('svmParam', 'var') && exist('mlrParam', 'var')))
    ap = zeros(param.numClass, length(param.C));
    mAP = 0;
  else
    try
      reuseComp = load(fileName);
    catch
      reuseComp.svmParam = [];
      reuseComp.modelSVM = [];
      reuseComp.ap = [];
    end

    for i = 1:param.numClass,
      %Prepare data here
      if(length(reuseComp.modelSVM) > i)
        modelSVM(i) = reuseComp.modelSVM(i);
        svmParam(i) = reuseComp.svmParam(i);
        ap(i, :) = reuseComp.ap(i, :);
      else
        trVal = 1;
        x = prepareDataSVM(X{trVal}, model.fNorm, model.iksvmN, model.isSparse, model.mapNanToZero);
        y = zeros(size(Y{trVal}));
        y(Y{trVal} == i) = 1; y(Y{trVal} ~= i) = -1;
        w = W{trVal};

        %Train all the SVM, in parallel for the cross validations
        tmpModel = trainSVM(x, y, w, param.C, param.w1, param.bias, param.weighted);
        % Train all the SVMs
        % for j = 1:length(param.C)
        %   tmpModel(j) = trainSVM(x, y, w, param.C(j), param.w1(j), param.bias(j), param.weighted);
        % end

        %Test all the SVMs on the validation data
        trVal = 2;
        x = prepareDataSVM(X{trVal}, model.fNorm, model.iksvmN, model.isSparse, model.mapNanToZero);
        y = zeros(size(Y{trVal}));
        y(Y{trVal} == i) = 1; y(Y{trVal} ~= i) = -1;
        w = W{trVal};
        for j = 1:length(param.C),
          [score pred] = testSVM(x, tmpModel(j));
          apI(j) = evalMetric(y, score, w, param.metric);
        end
        ap(i,:) = apI;
        fprintf('OVA AP for class %d.\n', i);
        ap(i,:)

        %Pick and store the best svm models
        [~, ind] = max(ap(i,:));
        modelSVM(i) = tmpModel(ind);
        svmParam(i) = ind;
      end

      %Compute the features for logistic regression
      trVal = 1;
      x = prepareDataSVM(X{trVal}, model.fNorm, model.iksvmN, model.isSparse, model.mapNanToZero);
      [lFeatures{trVal}(i,:), ~] = testSVM(x, modelSVM(i));
      
      trVal = 2;
      x = prepareDataSVM(X{trVal}, model.fNorm, model.iksvmN, model.isSparse, model.mapNanToZero);
      [lFeatures{trVal}(i,:), ~] = testSVM(x, modelSVM(i));
      
      %fileName = sprintf('%d-%s.mat', getPID(), datestr(now, 'yyyymmddTHHMMSS'));
      save(fileName, 'modelSVM', 'svmParam', 'ap');
    end
    fprintf('Trained OVA SVMs.\n\n');

    clear tmpModel
    keyboard;
    %Train a Multi class logistic on the features from the one versus all classifiers.
    for i = 1:length(param.CMLR),
      trVal = 2;
      tmpModel(i) = trainMCL(sparse(lFeatures{trVal}), Y{trVal}, W{trVal}, param.CMLR(i), param.biasMLR(i), param.mclWtType); 
      trVal = 1;
      [prob, pred] = testMCL(sparse(lFeatures{trVal}), tmpModel(i));
      mAP(i) = evalMetric(Y{trVal}, prob, W{trVal}, param.mclMetric);
    end
    mAP
    fprintf('Trained MCL.\n');
    [gr, ind] = max(mAP);
    mlrModel = tmpModel(ind);
    mlrParam = ind;
  end

  %Retrain on whole data or just the training data...
  clear score;
  if(param.useValForTraining)
    XX = [X{1}, X{2}];
    WW = [W{1}, W{2}];
    YY = [Y{1}, Y{2}];
  else
    XX = [X{1}];
    WW = [W{1}];
    YY = [Y{1}];
  end

  x = prepareDataSVM(XX, model.fNorm, model.iksvmN, model.isSparse, model.mapNanToZero);
  w = WW;

  if(param.weighted == false && ~issparse(x))
    yy = YY;
    y = zeros(param.numClass, numel(yy));
    for i = 1:param.numClass,
      y(i, yy == i) = 1;
      y(i, yy ~= i) = -1;
      Cs(i) = param.C(svmParam(i));
      w1s(i) = param.w1(svmParam(i));
      biass(i) = param.bias(svmParam(i));
    end
    svmModel = trainSVMV3(x, y, w, Cs, w1s, biass, param.weighted);
    fprintf('Trained OVA SVMs on whole data.\n');
    save(fileName, 'svmModel', 'score');
  
    wModel = cat(1, svmModel(:).w);
    biasModel = cat(1, svmModel(:).bias);
    assert(all(biasModel(1) == biasModel));
    if(biasModel(1) > 0)
      score = wModel(:,1:end-1)*x;
      score = bsxfun(@plus, score, wModel(:, end).*biasModel(1));
    else
      score = wModel(:,1:end)*x;
    end
    fprintf('Computed scores on the train set for training the MCL.\n');

  else
    for i = 1:param.numClass,
      yy = YY;
      y = zeros(size(yy));
      y(yy == i) = 1; y(yy ~= i) = -1;
      svmModel(i) = trainSVM(x, y, w, param.C(svmParam(i)), param.w1(svmParam(i)), param.bias(svmParam(i)), param.weighted);
      score(i,:) = testSVM(x, svmModel(i));
      save(fileName, 'svmModel', 'score');
    end
    fprintf('Trained OVA SVMs on whole data.\n');
  end

  %Train the Multiclass logistic here
  mlrModel = trainMCL(sparse(score), YY, w, param.CMLR(mlrParam), param.biasMLR(mlrParam), param.mclWtType); 
  fprintf('Trained MCL on whole data.\n');

  model.ap = ap;
  model.mAP = mAP;
  model.svmModel = svmModel;
  model.mlrModel = mlrModel;
  model.param = param;
  model.svmParamInd = svmParam;
  model.mlrParamInd = mlrParam;
end
