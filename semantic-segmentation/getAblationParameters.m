function PP = getAblationParameters(typ)
	typ
	%entryLevelNumClass = 40;
	
	switch typ,			
		case 'categorySpecific-level1-all',
		% Full system ablation 
			PP.featureTyp = {'gTexton', 'sift'};
			PP.selThresh = {[1], [1]};
			PP.featureCacheName = 'categorySpecific-all';
			PP.nVar = 1900;
			%PP.numClass = entryLevelNumClass;
		
		case 'categorySpecific-level1-gTexton',
		% Full system ablation 
			PP.featureTyp = {'gTexton'};
			PP.selThresh = {[1]};
			PP.featureCacheName = 'categorySpecific-gTexton';
			PP.nVar = 900;
			%PP.numClass = entryLevelNumClass;
		
		case 'categorySpecific-level1-sift',
		% Full system ablation 
			PP.featureTyp = {'colorSift'};
			PP.selThresh = {[1]};
			PP.featureCacheName = 'categorySpecific-sift';
			PP.nVar = 1000;
			%PP.numClass = entryLevelNumClass;

		case 'onlyApp',
		% Full system ablation 
			PP.featureTyp = {'generic', 'categorySpecific-sift'};
			PP.selThresh = {[1 2], [1 2]};
			PP.featureCacheName = 'ancOnlyApp';
			PP.selF = [42+[4 8 12 16 20], 68:83];
			PP.selF = [PP.selF, 101+PP.selF, 203:282];
			PP.nVar = length(PP.selF);
			%PP.numClass = entryLevelNumClass;

		case 'onlyGeom',
		% Full system ablation 
			PP.featureTyp = {'generic', 'categorySpecific-gTexton'};
			PP.selThresh = {[1 2], [1 2]};
			PP.featureCacheName = 'ancOnlyGeom';
			PP.selF = setdiff(1:101, [42+[4 8 12 16 20], 68:83]);
			PP.selF = [PP.selF, 101+PP.selF, 203:282];
			PP.nVar = length(PP.selF)
			%PP.numClass = entryLevelNumClass;

		case 'full',
		% Full system 
			PP.featureTyp = {'generic', 'categorySpecific-all'};
			PP.selThresh = {[1 2], [1 2]};
			PP.nVar = 282;
			PP.featureCacheName = 'ancFull';
			%PP.numClass = entryLevelNumClass;

		case 'sift+gTextons',
		% Only using the sift features 
			PP.featureTyp = {'colorSift', 'gTexton'};
			PP.selThresh = {[1 2], [1 2]};
			PP.nVar = 2*(1000+900);
			PP.featureCacheName = 'ancSiftGTexton';
			%PP.numClass = entryLevelNumClass;

		case 'gTextons',
		% Only using the sift features 
			PP.featureTyp = {'colorSift'};
			PP.selThresh = {[1 2]};
			PP.nVar = 2*900;
			PP.featureCacheName = 'ancGTexton';
			%PP.numClass = entryLevelNumClass;

		case 'sift',
		% Only using the sift features 
			PP.featureTyp = {'colorSift'};
			PP.selThresh = {[1 2]};
			PP.nVar = 2*1000;
			PP.featureCacheName = 'ancSift';
			%PP.numClass = entryLevelNumClass;
	
		case 'noAmodal',
			% No Amodal Completion
			PP.featureTyp = {'generic', 'categorySpecific-all'};
			PP.selThresh = {[1], [1]};
			PP.nVar = 282/2;
			PP.featureCacheName = 'noAncFull';
			%PP.numClass = entryLevelNumClass;

		case 'generic',
			% Only generic features
			PP.featureTyp = {'generic'};
			PP.selThresh = {[1 2]};
			PP.nVar = 202;
			PP.featureCacheName = 'ancOnlyGeneric';
			%PP.numClass = entryLevelNumClass;

		case 'categorySpecific',
			% Only generic features
			PP.featureTyp = {'categorySpecific-all'};
			PP.selThresh = {[1 2]};
			PP.nVar = 80;
			PP.featureCacheName = 'ancOnlyCategorySpecifc';
			%PP.numClass = entryLevelNumClass;

		case 'full+scene',
			% Full system 
			PP.featureTyp = {'generic', 'categorySpecific-all', 'scene'};
			PP.selThresh = {[1 2], [1 2], [1]};
			PP.nVar = 2*(101+40)+10;
			PP.featureCacheName = 'ancFullScene';
			%PP.numClass = entryLevelNumClass;


		case 'fullSC+scene',
			%Full + Scene
			PP.featureInDir = {amodalDir,amodalDir,amodalDir};
			PP.featureInExt = {'V5',v24, vScene};
			PP.nVar = 292;
			PP.selThresh = [1 2];
			PP.featureCacheName = 'anc';
			%PP.numClass = 4;
		
		case 'fullSC',
		% Full system 
			PP.featureInDir = {amodalDir,amodalDir};
			PP.featureInExt = {'V5',v24};
			PP.nVar = 282;
			PP.selThresh = [1 2];
			PP.featureCacheName = 'anc';
			%PP.numClass = 4;

		case 'full+sceneSVM',
			%Full + Scene
			PP.featureInDir = {amodalDir,amodalDir,amodalDir};
			PP.featureInExt = {'V5',v24,'V31'};
			PP.nVar = 293;
			PP.selThresh = [1 7];
			PP.featureCacheName = 'anc';
			%PP.numClass = entryLevelNumClass;
	end

end
