function x = prepareDataSVM(X, fNorm, iksvmN, sp, mapNanToZero)
% function x = prepareDataSVM(X, fNorm, iksvmN, sp, mapNanToZero)
%	sp is whether to return sparse or full output...
	if(mapNanToZero)
		X(isnan(X)) = 0;
	end
	x = normalizeFeatures(X, fNorm);
	if(iksvmN > 0)
		x = vl_homkermap(full(x),iksvmN);
	end
	if(sp)
		x = sparse(x);
	else
		x = full(x);
	end
end


