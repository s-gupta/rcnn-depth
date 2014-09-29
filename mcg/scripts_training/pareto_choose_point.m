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

% This function defines all parameters to tune in the training of MCG
% params = nyud_params(); %sf_mUCM_multi_3sc_u_4r_12k_params();

%% Compute all points from combining these base hierarchies
% The order of combination is hard-coded, see pareto_learning to change it
% The result will be loaded if already computed
function pareto_chose_point(params, num_regions)
  pareto = pareto_learning(params);

  %% Plot pareto front results

  % Play with the following two number to get the desired point

  % Which measure used to get the Pareto
  pareto_meas = 2;  % Jaccard at class level

  % Which particular point in the curve we choose
  pareto_id = -1;   % Around 14k candidates

  if(pareto_id == -1)
    pareto_id = find( pareto.st{pareto_meas}{params.n_r_cand}.mean_n_masks > num_regions, 1, 'first');
  end

  figure;
  x_limit = 1e6;
  titles{1} = {'Maximum achievable quality',['(Pascal segvoc ' params.gt_set_pareto ')'], 'Jaccard at instance level (Ji)'};
  titles{2} = {'Maximum achievable quality',['(Pascal segvoc ' params.gt_set_pareto ')'], 'Jaccard at class level (Jc)'};
  measures  = {'jaccard_object', 'jaccard_class'};
  colors = {'r','g','b','k','m','c'};
  for jj=1:length(measures)
      subplot(1,2,jj); hold on;

      % Singletons, pairs, and triplets
      % Pareto on the same measure that evaluated
      for ii=1:params.n_r_cand
          plot(pareto.st{jj}{ii}.mean_n_masks,pareto.st{jj}{ii}.(measures{jj}),[colors{ii} '-'])
      end
      
      for ii=1:params.n_r_cand
          % Pareto on one measure, evaluated on the other
          kk = setdiff([1 2],jj);
          plot(pareto.st{kk}{ii}.mean_n_masks,pareto.st{kk}{ii}.(measures{jj}),[colors{ii} '--'])
      end
      
      % Plot working point
      tmp = pareto.st{pareto_meas}{params.n_r_cand}.(measures{jj});
      plot(pareto.st{pareto_meas}{params.n_r_cand}.mean_n_masks(pareto_id),tmp(pareto_id),'k*')
      text(70,0.85,{[num2str(pareto.st{pareto_meas}{params.n_r_cand}.mean_n_masks(pareto_id)) ' c/i'],...
                    ['J=' num2str(tmp(pareto_id))]})
      
      % Make plot nicer
      title(strrep(titles{jj},'_','\_'))
      if (params.n_r_cand==3)
          legend({'Singletons same','Pairs same','Triplets same','Singletons diff','Pairs diff','Triplets diff'},4)
      elseif (params.n_r_cand==4)
          legend({'Singletons same','Pairs same','Triplets same','4-tuples','Singletons diff','Pairs diff','Triplets diff','4-tuples diff'},4)
      end
      grid minor
      grid on
      axis([10,x_limit,0.3,0.9])
      set(gca,'XScale','log')
  end

  %% Once you have found the desired point, save the parameters to file, i.e.,
  % the number of candidates from each hierarchy and for singletons, pairs, etc.
  sel_parameters  = pareto.pars{pareto_meas}{params.n_r_cand}(:,pareto_id);
  sampled_cands = pareto.base_stats{pareto_meas}.n_cands;
  total_cands = 0;
  clear n_cands;
  for ii=1:length(sel_parameters)
      if sel_parameters(ii)==0
          n_cands(ii) = 0; %#ok<*SAGROW>
      else
          n_cands(ii) = sampled_cands(sel_parameters(ii));
          total_cands = total_cands + pareto.base_stats{ii}.mean_n_masks(sel_parameters(ii));
      end
  end
  assert(abs(total_cands-pareto.st{pareto_meas}{params.n_r_cand}.mean_n_masks(pareto_id))<1e-5)
  n_cands = reshape(n_cands,length(params.hiers),params.n_r_cand)';
  n_cands

  save(params.files.pareto_point,'n_cands','pareto_id');
  disp(['Saved: ' params.files.pareto_point])
end
%% Write the curves to file (to put them in the paper)
% out_dir = '/Users/jpont/Publications/2014_CVPR/LaTeX/data/obj_cands/';
% 
% % Write up to singletons, pairs, and triplets (w.r.t pareto_meas)
% write_jaccard_to_file(pareto.st{pareto_meas}{1},fullfile(out_dir,[params.gt_set_pareto '_' params.hiers_id '_pareto_singletons.txt']));
% write_jaccard_to_file(pareto.st{pareto_meas}{2},fullfile(out_dir,[params.gt_set_pareto '_' params.hiers_id '_pareto_pairs.txt']));
% write_jaccard_to_file(pareto.st{pareto_meas}{3},fullfile(out_dir,[params.gt_set_pareto '_' params.hiers_id '_pareto_triplets.txt']));
% write_jaccard_to_file(pareto.st{pareto_meas}{4},fullfile(out_dir,[params.gt_set_pareto '_' params.hiers_id '_pareto_4tuples.txt']));
% 
% % Write selected pareto point
% % params = sf_mUCM_sub_params('sub1');
% load(params.files.pareto_point)
% fid = fopen(fullfile(out_dir,[params.gt_set_pareto '_' params.hiers_id '_pareto_selected_point.txt']),'w');
% fprintf(fid, 'ncands\tjac_class\tjac_instance\n');
% fprintf(fid, '%d\t%f\t%f\n', [pareto.st{pareto_meas}{4}.mean_n_masks(pareto_id); pareto.st{pareto_meas}{4}.jaccard_class(pareto_id); pareto.st{pareto_meas}{4}.jaccard_object(pareto_id)]);
% fclose(fid);
