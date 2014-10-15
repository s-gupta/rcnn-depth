function [feature sp spArea gt] = loadExamples(imList, paths, param, toTrain)
	gt = cell(1,length(imList));
	sp = cell(1,length(imList));
	feature = cell(1, length(imList));
	spArea = cell(1, length(imList));

	featureParam = param.featureParam;


	parfor i = 1:length(imList),
		fprintf('.');
		imName = imList{i};

		[iF iSP iSpArea] = getSPFeatures(imName, paths, featureParam);
		
		spArea{i}(1,:) = iSpArea;
		sp{i} = iSP;
		if(isfield(featureParam, 'selF'))
			feature{i} = iF(featureParam.selF,:);	
		else
			feature{i} = iF;
		end

		%Load the ground truth and project it to the superpixels
		if(toTrain)
			gtParam = param.gtParam;
			gtLabel = getGroundTruth(imName, 'segmentation', gtParam.classMapping);
			
			spHist = cell2mat(accumarray(iSP(:), gtLabel(:), [], @(x){linIt(histc(x,0:gtParam.numClass))})');
			spHist = bsxfun(@times, spHist, 1./(1+sum(spHist,1)));
			spHist = spHist(2:end,:);
			[maxVal, maxInd] = max(spHist, [], 1);
			gt{i} = maxInd;
			IND = maxVal >= gtParam.spFraction & spArea{i} >= gtParam.spArea;

			gt{i} = gt{i}(:, IND);
			feature{i} = feature{i}(:, IND);
			spArea{i} = spArea{i}(:, IND);
		end
	end
end
