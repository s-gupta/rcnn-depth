function model = trainMCL(x, y, w, C, bias, isWeighted) 
	if(isWeighted == false),
		model = train(double(y'), sparse(x), ...
			sprintf('-B %f -c %f -s %d -q',bias, C, 0),'col');
	elseif(isWeighted == true),
		model = trainWt(double(w'), double(y'), sparse(x), ...
			sprintf('-B %f -c %f -s %d -q',bias, C, 0),'col');
	end
	model.isSparse = true;
	model.isWeighted = isWeighted;
end
