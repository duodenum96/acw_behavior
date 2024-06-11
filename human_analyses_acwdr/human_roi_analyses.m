%% Electrode Specific Analyses
clear, clc, close all
cd C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\daviddata\
load("acwdata_acwdr.mat")
load('rest\001_rest_reordered.mat')
load('electrode_locations.mat')

el_locs = fields(elloc);
el_locs_names = {'Left Frontal', 'Right Frontal', 'Left Central', ...
    'Right Central', 'Left Parietal', 'Right Parietal', 'Left Temporal', ...
    'Right Temporal', 'Left Occipital', 'Right Occipital'};
%% Get some data
% For every subject, average rest runs, average task runs and take the
% difference
nchan = 64;
ntp = 47;

acws_rest_persubj = zeros(nsubj, nchan);
acws_self_persubj = zeros(nsubj, nchan);
acws_other_persubj = zeros(nsubj, nchan);

acws_rest_x_scalp = zeros(nsubj, ntp);
acws_self_x_scalp = zeros(nsubj, ntp);
acws_other_x_scalp = zeros(nsubj, ntp);

acws_sd_rest_persubj = zeros(nsubj, nchan);

acws_sd_rest_x_scalp = zeros(nsubj, ntp);

for i = 1:nsubj
    acws_rest_persubj(i, :) = mean(rest_acws.acwdrs_rest{i}, 2);
    acws_self_persubj(i, :) = mean(self_acws.acwdrs_self{i}, 2);
    acws_other_persubj(i, :) = mean(other_acws.acwdrs_other{i}, 2);

    acws_rest_x_scalp(i, :) = mean(rest_acws.acwdrs_rest{i}, 1);
    acws_self_x_scalp(i, :) = mean(self_acws.acwdrs_self{i}, 1);
    acws_other_x_scalp(i, :) = mean(other_acws.acwdrs_other{i}, 1);

    acws_sd_rest_persubj(i, :) = std(rest_acws.acwdrs_rest{i}, [], 2);

    acws_sd_rest_x_scalp(i, :) = std(rest_acws.acwdrs_rest{i}, [], 1);
end

restself_mean = ((acws_self_persubj - acws_rest_persubj) ./ acws_rest_persubj)*100;
restother_mean = ((acws_other_persubj - acws_rest_persubj) ./ acws_rest_persubj)*100;

restself_ave = mean(restself_mean);
restother_ave = mean(restother_mean);

rest_sd_ave = mean(acws_sd_rest_persubj);

%% Correct p values for multiple comparisons
all_pvals = zeros(3, length(el_locs));
for j = 1:length(el_locs)
    i_elloc = el_locs(j);
    chan_idx = elloc.(i_elloc{1});
    i_nchan = length(chan_idx);

    i_acws_rest_persubj = acws_rest_persubj(:, chan_idx);
    i_acws_self_persubj = acws_self_persubj(:, chan_idx);
    i_acws_other_persubj = acws_other_persubj(:, chan_idx);

    n_acws = length(i_acws_rest_persubj(:));
    data = [i_acws_rest_persubj(:); i_acws_self_persubj(:); i_acws_other_persubj(:)]; 
    groups = [ones(n_acws, 1); 2*ones(n_acws, 1); 3*ones(n_acws,1)];
   
    
    % 1-2
    [~, ~, stats] = ranksum(i_acws_rest_persubj(:), i_acws_self_persubj(:));
    [pval1,r,U,rsum,h] = ranksum_effect_size(i_acws_rest_persubj(:), i_acws_self_persubj(:));
    
    % 2-3
    [~, ~, stats] = ranksum(i_acws_other_persubj(:), i_acws_self_persubj(:));
    [pval2,r,U,rsum,h] = ranksum_effect_size(i_acws_other_persubj(:), i_acws_self_persubj(:));
    
    % 1-3
    [~, ~, stats] = ranksum(i_acws_rest_persubj(:), i_acws_other_persubj(:));
    [pval3,r,U,rsum,h] = ranksum_effect_size(i_acws_rest_persubj(:), i_acws_other_persubj(:));
    
    all_pvals(:, j) = [pval1, pval2, pval3];
end

correct_p_restself = bonf_holm(all_pvals(1, :));
correct_p_selfother = bonf_holm(all_pvals(2, :));
correct_p_restother = bonf_holm(all_pvals(3, :));

correct_p = [correct_p_restself; correct_p_selfother; correct_p_restother];
correct_p(correct_p > 1) = 1;

%% Rest  - Task Difference

for j = 1:length(el_locs)
    i_elloc = el_locs(j);
    chan_idx = elloc.(i_elloc{1});
    i_nchan = length(chan_idx);
    
    allcolors = slanCM('tab20', 21);
    col_1subj = repmat((1:21)', 1, i_nchan);
    col_all = [col_1subj(:); col_1subj(:); col_1subj(:)];
    colmat_all = zeros(length(col_all), 3);
    for i = 1:length(col_all)
        i_subj = col_all(i);
        colmat_all(i, :) = allcolors(i_subj, :);
    end

    i_acws_rest_persubj = acws_rest_persubj(:, chan_idx);
    i_acws_self_persubj = acws_self_persubj(:, chan_idx);
    i_acws_other_persubj = acws_other_persubj(:, chan_idx);

    n_acws = length(i_acws_rest_persubj(:));
    data = [i_acws_rest_persubj(:); i_acws_self_persubj(:); i_acws_other_persubj(:)]; 
    groups = [ones(n_acws, 1); 2*ones(n_acws, 1); 3*ones(n_acws,1)];
    
    swarmchart(groups, data, [], colmat_all, 'filled')
    line([0.75 1.25], [median(data(groups==1)) median(data(groups==1))], ...
        'color', 'k', 'LineWidth', 1.5)
    line([1.75 2.25], [median(data(groups==2)) median(data(groups==2))], ...
        'color', 'k', 'LineWidth', 1.5)
    line([2.75 3.25], [median(data(groups==3)) median(data(groups==3))], ...
        'color', 'k', 'LineWidth', 1.5)
    set(gca, 'YGrid', 'on')
    
    % 1-2
    [~, ~, stats] = ranksum(i_acws_rest_persubj(:), i_acws_self_persubj(:));
    [~,r,U,rsum,h] = ranksum_effect_size(i_acws_rest_persubj(:), i_acws_self_persubj(:));
    stattext = get_stattext_wilcox(stats.zval, correct_p(1, j), r, 'horizontal');
    sigstar_duodenal({[1 2]}, correct_p(1, j),stattext, 0)
    % 2-3
    [~, ~, stats] = ranksum(i_acws_other_persubj(:), i_acws_self_persubj(:));
    [~,r,U,rsum,h] = ranksum_effect_size(i_acws_other_persubj(:), i_acws_self_persubj(:));
    stattext = get_stattext_wilcox(stats.zval, correct_p(2, j), r, 'horizontal');
    sigstar_duodenal({[2 3]}, correct_p(2, j),stattext, 0)
    % 1-3
    [~, ~, stats] = ranksum(i_acws_rest_persubj(:), i_acws_other_persubj(:));
    [~,r,U,rsum,h] = ranksum_effect_size(i_acws_rest_persubj(:), i_acws_other_persubj(:));
    stattext = get_stattext_wilcox(stats.zval, correct_p(3, j), r, 'horizontal');
    sigstar_duodenal({[1 3]}, correct_p(3, j),stattext, 0)
    ax = gca;
    % set(ax.Children(1), 'Position', [2.0414746s5437788,0.691159842021277,0])
    % set(ax.Children(3), 'Position', [3.004608294930875,0.60864295,0])
    % set(ax.Children(5), 'Position', [1.43778801843318,0.570773235372341,0])
    
    % kruskalwallis(da  ta, groups)
    
    xticks([1 2 3])
    xticklabels(["Rest", "Self", "Other"])
    ylabel("Mean of timescale (\tau) across windows")
    
    fontsize(gcf, 14, "points")
    title(el_locs_names(j))
    pause % adjust text manually
    saveas(gcf, sprintf('figs_acwdr\\channels\\acwdiff\\acwdiff_stat_acwdr_specific_topo_%s.fig', i_elloc{1}))
    saveas(gcf, sprintf('figs_acwdr\\channels\\acwdiff\\acwdiff_stat_acwdr_specific_topo_%s.jpg', i_elloc{1}))
    close all
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

