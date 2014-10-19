function param = getClassifierParam(typ, data)

switch(typ)
  case 'svm-categorySpecific',
    C = 10.^[-3:1:3, -3:1:3, -3:1:3];
    b = C; b(:) = 1;
    w = C; w(:) = 10; w(1:2*length(C)/3) = 3; w(1:length(C)/3) = 1;
    param = ...
      struct('iksvmN', 3, 'C', C, 'w1', w, 'bias', b, ...
        'numClass', 1, 'CMLR', 1, 'biasMLR', 1, ...
        'mapNanToZero', false, 'useValForTraining', data.useVal, 'weighted', false, ...
        'metric', 'ap1min', 'mclMetric', 'multiClassAccuracy', 'mclWeighted', true);

  case 'svm-full',
    param.C = 10.^[-2:1:1 -2:1:1];
    param.w1 = [1 1 1 1 3 3 3 3];
    param.bias = ones(size(param.w1));
    param.CMLR = 1;
    param.biasMLR = 1;
    param.iksvmN = 3; %iksvmN;
    param.numClass = data.numClass;
    param.metric = 'ap1';
    param.mapNanToZero = true;
    param.weighted = false; 
    param.mclMetric = 'multiClassAccuracy'; 
    param.mclWeighted = 1;
    param.useValForTraining = data.useVal;

  case 'tree-full',
    param.nTree = 40; %nTree;
    param.minLeaf = 20; %minLeaf;
    param.mapNaNToZero = true; %mapNaNToZero;
    param.useWeights = true; %useWeights;
    param.useValForTraining = data.useVal;
    param.nVarToSample = round(data.nVar/3);
    param.seed = 0;
    param.useVal = data.useVal;

  case 'svm-scene',
    param.C = 10.^[-6:1:1 -6:1:1]; % -3:1:1];
    param.w1 = [1 1 1 1 1 1 1 1, 3 3 3 3 3 3 3 3]; %, 10 10 10 10 10 10];
    param.bias = ones(size(param.w1));
    param.CMLR = [1];
    param.biasMLR = [1];
    param.iksvmN = 3;
    param.numClass = data.numClass;
    param.metric = 'ap1';
    param.mclMetric = 'multiClassAccuracy';
    param.mapNanToZero = true;
    param.weighted = false; 
    param.useValForTraining = data.useVal;
end



end
