args = {};
if strcmp(jobName, 'compute_ss_features')
  imset = {'train', 'val', 'test'};
  for i = 1:length(imset),
    imlist = getImageSet(imset{i});
    for j = 1:length(imlist),
      args{end+1} = {imlist{j}};
    end
  end
  jobParam = struct('numThreads', 1, 'codeDir', pwd(), 'preamble', '', 'matlabpoolN', 0, 'globalVars', {{}}, 'fHandle', @script_ss, 'numOutputs', 0);
  resourceParam = struct('mem', 8, 'hh', 20, 'numJobs', 34, 'ppn', 1, 'nodes', 1, 'logDir', '/work4/sgupta/pbsBatchDir/', 'queue', 'psi', 'notif', false, 'username', 'sgupta', 'headNode', 'zen');
  keyboard;
  [jobId jobDir] = jobParallel(jobName, resourceParam, jobParam, args);
end

if strcmp(jobName, 'detectors_for_sseg')
  res = rcnn_all('task-detection', 'rgb_hha', 'trainval', 'test');
  res = rcnn_all('task-detection', 'rgb_hha', 'train', 'val');
  res = rcnn_all('task-detection', 'rgb_hha', 'val', 'train');
end

if strcmp(jobName, 'copy_boxes')
  p = get_paths();
  
  s = sprintf('cp %s/detections/*boxes*test*.mat %s/.', fullfile(p.detection_dir, 'detector', 'rgb_hha_30000_trainval'), p.ss_det_dir); 
  fprintf('%s\n\n', s); system(s);
  s = sprintf('cp %s/pr-curves/*pr*test*.mat %s/.', fullfile(p.detection_dir, 'detector', 'rgb_hha_30000_trainval'), p.ss_det_dir); 
  fprintf('%s\n\n', s); system(s);
  
  s = sprintf('cp %s/detections/*boxes*train*.mat %s/.', fullfile(p.detection_dir, 'detector', 'rgb_hha_30000_val'), p.ss_det_dir); 
  fprintf('%s\n\n', s); system(s);
  s = sprintf('cp %s/pr-curves/*pr*train*.mat %s/.', fullfile(p.detection_dir, 'detector', 'rgb_hha_30000_val'), p.ss_det_dir); 
  fprintf('%s\n\n', s); system(s);
  
  s = sprintf('cp %s/detections/*boxes*val*.mat %s/.', fullfile(p.detection_dir, 'detector', 'rgb_hha_30000_train'), p.ss_det_dir); 
  fprintf('%s\n\n', s); system(s);
  s = sprintf('cp %s/pr-curves/*pr*val*.mat %s/.', fullfile(p.detection_dir, 'detector', 'rgb_hha_30000_train'), p.ss_det_dir); 
  fprintf('%s\n\n', s); system(s);
end

if strcmp(jobName, 'deepDetFeatures')
  args = {};
  imSet = {'train', 'val', 'test'};
  typ = 'detection-box';
  catName = getMetadata('task-detection'); catName = catName.name;
  detDir = p.ss_det_dir; nmsThresh = 0.3;
  assignTyp = {'bestScore'}; pThresh = 0.50;
  threshSet = 'val'; boxFieldName = 'boxes'; out_dir = p.ss_det_feature;
  K = 1; 
  for i = 2:length(imSet), for j = 1:length(assignTyp), for k = 1:length(pThresh), for l = 1:K,
    imlist = getImageSet(imSet{i});
    args{end+1} = {imSet{i}, typ, l:K:length(imlist), catName, detDir, assignTyp{j}, pThresh(k), threshSet, nmsThresh, boxFieldName, p.ss_generic, out_dir};
    compute_detection_features(args{end}{:});
  end, end, end, end

  % jobParam = struct('numThreads', 1, 'codeDir', pwd(), 'preamble', '', 'matlabpoolN', 4, 'globalVars', {{}}, 'fHandle', @compute_detection_features, 'numOutputs', 0);
  % resourceParam = struct('mem', 4, 'hh', 8, 'numJobs', -1, 'ppn', 4, 'nodes', 1, 'logDir', '/work4/sgupta/pbsBatchDir/', 'queue', 'psi', 'notif', false, 'username', 'sgupta', 'headNode', 'zen');
  % [jobId jobDir] = jobParallel(jobName, resourceParam, jobParam, args);
end

if strcmp(jobName, 'ablationStudy')
	classMapping = 'classMapping40';
	trSet = 'train1'; valSet = 'train2';	testSet = 'val';
	classifier = {'svm-full'};
	typ = {'full+deepdet', 'full', 'generic', 'categorySpecific'};

	for i = 1:length(classifier),
		for j = 1, %:length(typ), 
			evalResAblation(i,j) = wrapperTrainTestBenchmarkModel(trSet, valSet, testSet, typ{j}, classifier{i}, classMapping);
			fprintf('%s - %s, accuracy found to be, %0.3f.\n', classifier{i}, typ{j}, evalResAblation(i,j).fwavacc);
		end
	end
end


if strcmp(jobName, 'full-sseg')
	classMapping = 'classMapping40';
	trSet = 'train'; valSet = 'val';	testSet = 'test';
	classifier = {'svm-full'};
	typ = {'full', 'full+deepdet'};

	for i = 1:length(classifier),
		for j = 1:2, %length(typ), 
			evalRes(i,j) = wrapperTrainTestBenchmarkModel(trSet, valSet, testSet, typ{j}, classifier{i}, classMapping);
			fprintf('%s - %s, accuracy found to be, %0.3f.\n', classifier{i}, typ{j}, evalRes(i,j).fwavacc);
		end
	end
end


