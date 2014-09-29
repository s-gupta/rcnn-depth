% ------------------------------------------------------------------------ 
%  Copyright (C)
%  Universitat Politecnica de Catalunya BarcelonaTech (UPC) - Spain
% 
%  Jordi Pont-Tuset <jordi.pont@upc.edu>
%  June 2013
% ------------------------------------------------------------------------ 

function [pareto_stats, base_params] = pareto_combination(base_stats, comb_sequence, measure)

n_base_stats = length(base_stats);
n_merges = size(comb_sequence,1);

assert(n_base_stats==n_merges+1);
assert(size(comb_sequence,2)==3);

params = cell(n_base_stats+n_merges,1);
pareto_stats = cell(n_base_stats+n_merges,1);

for ii=1:n_base_stats
    pareto_stats{ii} = base_stats{ii};
end
for ii=1:size(comb_sequence,1)
    disp([ii size(comb_sequence,1)]);
    [~,pareto_stats{comb_sequence(ii,3)}, params{comb_sequence(ii,3)}] =...
          combine_masks(pareto_stats{comb_sequence(ii,1)},pareto_stats{comb_sequence(ii,2)}, measure);
end

% Allocate
base_params = cell(n_base_stats+n_merges,1);

% Populate leave parameterizations
for ii=1:n_base_stats
    n_results = size(base_stats{ii}.max_J,2);
    base_params{ii} = zeros(n_base_stats,n_results);
    for kk=1:n_results
        base_params{ii}(ii,kk) = kk;
    end
end

% Evolve through combination sequence
for ii=1:n_merges
    assert(comb_sequence(ii,3)==n_base_stats+ii)
    n_results = length(params{n_base_stats+ii});
    base_params{n_base_stats+ii} = zeros(n_base_stats,n_results);
    for kk=1:n_results
        merged1 = params{n_base_stats+ii}(kk).params1;
        merged2 = params{n_base_stats+ii}(kk).params2;
        if merged1<size(base_params{comb_sequence(ii,1)},2)+1
            pars1 = base_params{comb_sequence(ii,1)}(:,merged1);
        elseif merged1==size(base_params{comb_sequence(ii,1)},2)+1
            pars1 = zeros(n_base_stats,1);
        else
            error('Oh oh')
        end
        if merged2<size(base_params{comb_sequence(ii,2)},2)+1
            pars2 = base_params{comb_sequence(ii,2)}(:,merged2);
        elseif merged2==size(base_params{comb_sequence(ii,2)},2)+1
            pars2 = zeros(n_base_stats,1);
        else
            error('Oh oh')
        end
        if sum(pars1.*pars2)>0
            error('Oh oh 2')
        end
        base_params{n_base_stats+ii}(:,kk) = pars1 + pars2;
    end
end

