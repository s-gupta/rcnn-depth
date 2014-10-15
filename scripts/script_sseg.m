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
