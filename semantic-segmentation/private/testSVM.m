function [score, pred] = testSVM(x, model);
	isSparse = model.isSparse;
	isWeighted = model.isWeighted; 

	if(isSparse == true)
		[pred, ~, score] = predict(zeros(1,size(x,2))', x, ...
			model, '', 'col');
	elseif(isSparse == false)
		[pred, ~, score] = predictDense(zeros(1,size(x,2))', x, ...
			model, '', 'col');
	end
	pred = pred';
	score = score';
end
