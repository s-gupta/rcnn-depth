function model = trainSVM(x, y, w, C, w1, bias, isWeighted)
  isSparse = issparse(x);

  %% Hack used only for parallel training...
  for i = 1:length(C),
    args{i} = sprintf('-s 3 -B %f -c %f -w1 %f', ...
        bias(i), C(i), w1(i));
  end

  if(isSparse == true & isWeighted == true)
    model = trainWt(double(w')./max(w'), double(y'), x, ... 
      sprintf('-s 3 -B %f -c %f -w1 %f', ...
        bias, C, w1),...
        'col');
    model = fixFlipping(model);

  elseif(isSparse == true & isWeighted == false)
    model = train(double(y'), x,...
      sprintf('-s 3 -B %f -c %f -w1 %f', ...
        bias, C, w1),...
        'col');
    model = fixFlipping(model);

  elseif(isSparse == false & isWeighted == true)
    model = trainWt(double(w')./max(w'), double(y'), sparse(x), ... 
      sprintf('-s 3 -B %f -c %f -w1 %f', ...
        bias, C, w1),...
        'col');
    model = fixFlipping(model);
    fprintf('Unimplemented!!');

  elseif(isSparse == false & isWeighted == false)
    if(length(args) == 1)
      model = trainDense(double(y'), x, ... 
        sprintf('-s 3 -B %f -c %f -w1 %f', ...
          bias, C, w1),...
          'col');
      model = fixFlipping(model);
    else
      fprintf('Using openMP SVM training.\n');
      if(size(y,1) == 1)
        model = trainDenseOpenMP_400(double(y'), x, args, 'col');
      else
        model = trainDenseOpenMPY(double(y'), x, args, 'col');
      end
      for i = 1:length(model),
        m(i) = fixFlipping(model{i}); 
      end
      model = m;
    end
  end

  for i = 1:length(model),
    model(i).isSparse = isSparse;
    model(i).isWeighted = isWeighted;
  end
end

function model = fixFlipping(model)
  model.w = model.w.*model.Label(1);
  model.Label = model.Label.*model.Label(1);
end
