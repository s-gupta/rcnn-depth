% ------------------------------------------------------------------------ 
%  Copyright (C)
%  Universitat Politecnica de Catalunya BarcelonaTech (UPC) - Spain
%  University of California Berkeley (UCB) - USA
% 
%  Jordi Pont-Tuset <jordi.pont@upc.edu>
%  Pablo Arbelaez <arbelaez@berkeley.edu>
%  June 2014
% ------------------------------------------------------------------------ 
% This file is part of the MCG package presented in:
%    Arbelaez P, Pont-Tuset J, Barron J, Marques F, Malik J,
%    "Multiscale Combinatorial Grouping,"
%    Computer Vision and Pattern Recognition (CVPR) 2014.
% Please consider citing the paper if you use this code.
% ------------------------------------------------------------------------
function out = cache_mcg_features(params, ucms, imName, cache_dir)
  if(isempty(ucms)), 
    for ii = 1:length(params.hier_dirs),
      tmp = load(fullfile_ext(params.hier_dirs{ii}, imName, 'mat')); 
      ucms{ii} = tmp.ucm2; 
    end
  end

  if(~isempty(cache_dir)), cache = true; cache_file = fullfile_ext(cache_dir, imName, 'mat');
  else cache = false; end
  
  n_hiers = length(params.hier_dirs);

  % Load Pareto parameters
  n_cands = loadvar(params.files.pareto_point,'n_cands');

  % Number of regions per candidate
  assert(n_hiers==size(n_cands,2));

  % Is it already computed?
  try
    assert(cache);
    out = load(cache_file);
  catch
    % Read all hierarchies
    lps = [];
    ms  = cell(n_hiers,1);
    ths = cell(n_hiers,1);
    all_ucms = [];
    for ii=1:n_hiers
        ucm2 = ucms{ii};
        
        all_ucms = cat(3,all_ucms, ucm2);

        % Read the UCM as a hierarchy
        curr_hier = ucm2hier(ucm2);
        ths{ii}.start_ths = curr_hier.start_ths';
        ths{ii}.end_ths = curr_hier.end_ths';
        ms{ii} = curr_hier.ms_matrix;
        lps = cat(3, lps, curr_hier.leaves_part);
    end

    % Get full cands, represented on a fused hierarchy
    tt = tic();
    [f_lp,f_ms,cands,start_ths,end_ths] = full_cands_from_hiers(lps,ms,ths,n_cands);
    fprintf('Time for full_cands_from_hiers: %0.3f\n', toc(tt));

    % Hole filling and complementary candidates
    tt = tic();
    [cands_hf, cands_comp] = hole_filling(double(f_lp), double(f_ms), cands); %#ok<NASGU>
    fprintf('Time for hole_filling: %0.3f\n', toc(tt));
    
    % Select which candidates to keep (Uncomment just one line)
    cands = cands_hf;                       % Just the candidates with holes filled
    % cands = [cands_hf; cands_comp];         % Holes filled and the complementary
    % cands = [cands; cands_hf; cands_comp];  % All of them
    
    % Compute base features
    tt = tic();
    b_feats = compute_base_features(f_lp, f_ms, all_ucms);
    fprintf('Time for compute_base_features: %0.3f\n', toc(tt));
    b_feats.start_ths = start_ths;
    b_feats.end_ths   = end_ths;
    b_feats.im_size   = size(f_lp);

    % Filter by overlap
    tt = tic();
    red_cands = mex_fast_reduction(cands-1,b_feats.areas,b_feats.intersections,params.J_th);
    fprintf('Time for mex_fast_reduction: %0.3f\n', toc(tt));

    % Compute full features on reduced cands
    tt = tic();
    [feats, bboxes] = compute_full_features(red_cands,b_feats);
    fprintf('Time for compute_full_features: %0.3f\n', toc(tt));

    tt = tic();
    candidates.superpixels = f_lp;
    candidates.labels = cands2labels(red_cands, f_ms);
    fprintf('Time for cand2labels: %0.3f\n', toc(tt));
    
    b_feats_intersection = sparse(b_feats.intersections);
    f_ms = sparse(f_ms);

    % Save
    % save(cache_file, 'f_ms', 'f_lp', 'candidates', 'feats', 'bboxes', 'red_cands', 'b_feats_intersection');
    out = struct('f_ms', f_ms, 'f_lp', f_lp, 'feats', feats, 'red_cands', red_cands, 'b_feats_intersection', b_feats_intersection, 'bboxes', bboxes);
    if(cache), save(cache_file, '-STRUCT', 'out'); end
  end
end
