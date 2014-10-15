function modelOut = trainTreeBagging(X, y, W, param)
 	% X is DxN
	% y is Nx1
	% W is Nx1
	if(param.useVal)
		X = cat(2, X{:});
		y = cat(2, y{:});
		W = cat(2, W{:});
	else
		X = X{1};
		y = y{1};
		W = W{1};
	end

	switch param.useWeights,
		case 0,
			W(:) = 1;
		case 1,
			W = W;
		case 2, 
			W = log(W+1);
	end
	if(param.mapNaNToZero)
		X(isnan(X)) = 0;
	end
	if(~isfield(param, 'OOBVarImp'))
		param.OOBVarImp = 'on';
	end

	options = statset('UseParallel','always', 'Streams', RandStream('mlfg6331_64','seed',param.seed),'UseSubStreams','always');
	%options = statset('UseParallel','never', 'Streams', RandStream('mlfg6331_64','seed',param.seed),'UseSubStreams','always');

	model = TreeBagger(param.nTree, X', y,'oobpred','on','weights',W, 'method','classification','priorprob','empirical','minleaf', param.minLeaf,'NPrint', 1,...
		'NVarToSample',param.nVarToSample,'OOBVarImp',param.OOBVarImp,'Options', options);

	for i = 1:length(model.ClassNames),
		order(i) = str2num(model.ClassNames{i});
	end
	modelOut.model = model;
	modelOut.order = order;
	modelOut.mapNaNToZero = param.mapNaNToZero;
end
