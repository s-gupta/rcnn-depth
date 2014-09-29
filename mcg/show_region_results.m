% Performance for our object categories
function show_region_results(result_dir, method, legends),
% Performance for our object categories

titles = {'mean(maxJ)', 'mean(maxJ > 0.5)', 'mean(maxJ > 0.75)'};

task = {'35 Object Classes from Gupta et al.', '21 Object Classes from Lin et al.'};

x_limit = 1e4;
cols = lines(length(method));
catg{1} = [3 4 5 6 7 9 10 11 12 13 14 15 16 17 18 19 20 21 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 39 40];

%%
i = 1;
  for j = 1:length(method),
    dt = load(sprintf('%s/%s-eval.mat', result_dir, method{j}));
    stats = dt.stats;
    n_cands = stats.mean_n_masks;
    meanmax = []; recall50 = []; recall75 = [];
    for k = 1:length(catg{i}),
      max_J = stats.max_J(stats.obj_classes == catg{i}(k), :);
      meanmax(k,:) = mean(max_J, 1);
      recall50(k,:) = mean(max_J > 0.5, 1);
      recall75(k,:) = mean(max_J > 0.75, 1);
    end
    meanmax = mean(meanmax);
    recall50 = mean(recall50);
    recall75 = mean(recall75);
    if(j <= 0),
      n_cands = max(n_cands);
      meanmax = max(meanmax);
      recall50 = max(recall50);
      recall75 = max(recall75);
    end
    a(j) = struct('meanmax', meanmax, 'recall50', recall50, 'recall75', recall75, 'n_cands', n_cands); 
  end
  % Plot these things!
  figure(i); clf;
  subplot(2,2,1); hold on;
  for j = 1:length(method),
    if(length(a(j).n_cands) == 1)
      plot(a(j).n_cands, a(j).meanmax, ['.'], 'Color', cols(j,:), 'MarkerSize', 20);
    else
      plot(a(j).n_cands, a(j).meanmax, ['-'], 'Color', cols(j,:), 'LineWidth', 1);
    end
  end
  
  subplot(2,2,2); hold on;
  for j = 1:length(method),
    if(length(a(j).n_cands) == 1)
      plot(a(j).n_cands, a(j).recall50, ['.'], 'Color', cols(j,:), 'MarkerSize', 20);
    else
      plot(a(j).n_cands, a(j).recall50, ['-'], 'Color', cols(j,:), 'LineWidth', 1);
    end
  end
  
  subplot(2,2,3); hold on;
  for j = 1:length(method),
    if(length(a(j).n_cands) == 1)
      plot(a(j).n_cands, a(j).recall75, ['.'], 'Color', cols(j,:), 'MarkerSize', 20);
    else
      plot(a(j).n_cands, a(j).recall75, ['-'], 'Color', cols(j,:), 'LineWidth', 1);
    end
  end
  
  for j = 1:3,
    subplot(2,2,j);
    title(sprintf('%s %s', task{i}, strrep(titles{j},'_','\_')));
    legend(legends, 4);
    grid minor
    grid on
    axis([10, x_limit, 0.25, 0.9])
    set(gca,'XScale','log')
  end
  
  %% Make the plots that we want to put in the paper
  figure(2+i); clf; hold on;
  for j = 1:length(method),
    plot(a(j).n_cands, a(j).meanmax, ['-'], 'Color', cols(j,:), 'LineWidth', 2.0);
    
  end
  set(gca, 'FontName', 'times');
  title(task{i}, 'FontSize', 16);
  h_legend = legend(legends, 4);
  set(h_legend,'FontSize', 16, 'FontName', 'times');
  grid minor
  grid on
  xlabel('Number of candidates', 'FontSize', 18);
  ylabel('Coverage (Average Jaccard Index over Classes)', 'FontSize', 18);
  axis([10, x_limit, 0.25, 0.9])
  set(gca,'XScale','log');
  set(gca, 'FontSize', 18);
  
  % saveas(gcf, sprintf('paper_figures/newnewbaseline-region-recall-%d.fig', i))
  % set(gcf,'paperposition',get(gcf,'position')/100);
  % saveas(gcf, sprintf('paper_figures/newnewbaseline-region-recall-%d.eps', i), 'epsc');
end
