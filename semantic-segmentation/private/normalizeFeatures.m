function [features, fNorm] = normalizeFeatures(features, fNorm)
	if(~exist('fNorm','var'))
		fNorm.minimum = nanmin(features')';	
		features = bsxfun(@minus,features,fNorm.minimum);
		fNorm.maximum = nanmax(features')';
		features = bsxfun(@times,features,1./(fNorm.maximum+eps));
	else
		features = bsxfun(@minus,features,fNorm.minimum);
		features = bsxfun(@times,features,1./(fNorm.maximum+eps));
	end
end
