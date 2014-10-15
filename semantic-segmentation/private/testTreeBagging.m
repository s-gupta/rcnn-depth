function prob = testTreeBagging(model, X)
 	% X is DxN
	if(isfield(model, 'mapNaNToZero') & model.mapNaNToZero)
		X(isnan(X)) = 0;
	end
	[A B] = model.model.predict(X');
	b = B';
	clear B
	B(model.order,:) = b;
	prob = B;
end
