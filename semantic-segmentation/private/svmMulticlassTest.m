function [score, prob, pred] = svmMulticlassTest(model, X)
%function [score, prob, pred] = svmMulticlassTest(model, X)
% Input: 
	x = prepareDataSVM(X, model.fNorm, model.iksvmN, model.isSparse, model.mapNanToZero);
	for i = 1:model.numClass,
		[score(i,:), gr] = testSVM(x, model.svmModel(i));
	end
	%Test the Multiclass logistic here
	[prob, pred] = testMCL(sparse(score), model.mlrModel);
end
