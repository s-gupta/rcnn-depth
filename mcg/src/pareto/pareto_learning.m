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

function pareto = pareto_learning(params)
n_hiers = length(params.hiers);

% Are the results already computed?
if exist(params.files.pareto_full, 'file')
    load(params.files.pareto_full)
    disp(['Loaded: ' params.files.pareto_full '.'])
    recompute = 0;
else
    disp(['RECOMPUTING: ' params.files.pareto_full '.'])
    recompute = 1;
end

if recompute
    % Store params.n_r_cand to output 
    pareto.n_r_cand = params.n_r_cand;
    
    % Get base stats: each of the hierarchies independently 
    % and for singletons, pairs, etc. #ok<*AGROW>
    for ii=1:n_hiers
        curr_stats = get_single_hier_stats(params,ii);
        for jj=1:params.n_r_cand
            pareto.base_stats{ii+n_hiers*(jj-1)} = curr_stats(jj);
        end
    end

    % Store where the base candidates are for each type of number of regions
    for jj=1:params.n_r_cand
        pareto.base_ids{jj} = n_hiers*(jj-1)+(1:n_hiers);
    end
    
    % Combination sequence: first all singletons, pairs, and triplets, etc.
    curr_n = n_hiers*params.n_r_cand;
    pareto.comb_sequence = [];
    if n_hiers==1
        pareto.global_ids = 1:params.n_r_cand;
    else  
        for jj=1:params.n_r_cand
            curr_n = curr_n + 1;
            pareto.comb_sequence = [pareto.comb_sequence;...
                                    pareto.base_ids{jj}(2) pareto.base_ids{jj}(1) curr_n];
            for ii=3:length(pareto.base_ids{jj})
                pareto.comb_sequence = [pareto.comb_sequence;...
                                        pareto.base_ids{jj}(ii) curr_n curr_n+1];
                curr_n = curr_n + 1;
            end

            % Store the id of the global combination of singletons, pairs, etc.
            pareto.global_ids(jj) = curr_n;
        end
    end
    
    % After that singletons-pairs, the result with triplets, etc.
    curr_n = curr_n + 1;
    pareto.comb_sequence = [pareto.comb_sequence;...
                            pareto.global_ids(1) pareto.global_ids(2) curr_n];
    for jj=3:params.n_r_cand
        pareto.comb_sequence = [pareto.comb_sequence;...
                                pareto.global_ids(jj) curr_n curr_n+1];
        curr_n = curr_n + 1;
    end

    % Compute pareto front
    pareto.measures = {'jaccard_object', 'jaccard_class'};
    parfor jj=1:length(pareto.measures)
        [pareto_stats{jj}, pars{jj}] = pareto_combination(pareto.base_stats, pareto.comb_sequence, pareto.measures{jj}); %#ok<AGROW>
    end

    % Save full Pareto combination stats
    % save(fullfile(res_dir,[params.gt_set_pareto '_pareto_stats_' res_id '.mat']), 'pareto_stats', 'params','base_stats', 'comb_sequence', 'id','times');

    % Save just the results at the Pareto front for
    % singletons, pairs, triplets, etc.
    for ii=1:length(pareto.measures)
        pareto.st{ii}{1}      = pareto_stats{ii}{pareto.global_ids(1)};
        pareto.pars{ii}{1}    = pars{ii}{pareto.global_ids(1),1};  % Just singletons
        for jj=2:params.n_r_cand
            pareto.st{ii}{jj}      = pareto_stats{ii}{end-params.n_r_cand+jj}; % Singletons + pairs + ...
            pareto.pars{ii}{jj}    = pars{ii}{end-params.n_r_cand+jj,1};
        end
    end
    if ~exist(fullfile(root_dir,'datasets',params.database,'pareto'),'dir')
        mkdir(fullfile(root_dir,'datasets',params.database,'pareto'))
    end
    
    save(params.files.pareto_full, 'pareto');
end
