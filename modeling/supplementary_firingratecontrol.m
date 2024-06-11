%% Analyze results from clustered network
clear, clc, close all
cd C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\modeling_3\
addpath C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\
addpath C:\Users\duodenum\Desktop\brain_stuff\misc\slanCM\

nsims=30;
colors = repmat(1:nsims, 360, 1);
colors = [colors(:); colors(:)];
colmat = zeros(length(colors), 3);
all_cols = slanCM('glasbey', 30);

for i = 1:length(colors)
    i_col = colors(i);
    colmat(i, :) = all_cols(i_col, :);
end

load cellreports_review\dtinetwork_control_fr.mat
% data dims: roi, window, beta, simulations

var_rest = squeeze(std(acwdrs_rest, [], 2));
var_task = squeeze(std(acwdrs_task, [], 2));

mean_rest = squeeze(mean(acwdrs_rest, 2));
mean_task = squeeze(mean(acwdrs_task, 2));

rt = ((mean_task - mean_rest) ./ mean_rest) * 100;
nshuffles = 7;

ps = zeros(1, nshuffles);
rs = zeros(1, nshuffles);
zs = zeros(1, nshuffles);
for i = 1:nshuffles
    i_mean_rest = mean_rest(:, i, :);
    i_mean_task = mean_task(:, i, :);
    [ps(i), ~, stats] = ranksum(i_mean_rest(:), i_mean_task(:));
    zs(i) = stats.zval;
    [~, rs(i)] = ranksum_effect_size(i_mean_rest(:), i_mean_task(:));
end
ps_correct = bonf_holm(ps, 0.05);
%% Mean Change
surrogates = [0 0.5 1 1.5 2 2.5 3];
shufflestring = string(surrogates);
close all
figure('Position',[-1.4606e+03 129.8000 784.8000 655.2000]);
tiledlayout(3,3, 'TileSpacing','tight', 'Padding','tight')
% tiledlayout(2, 2)
for i = 1:nshuffles
    i_mean_rest = mean_rest(:, i, :);
    i_mean_task = mean_task(:, i, :);

    i_rt = rt(:, i, :);
    i_var_rest = var_rest(:, i, :);

    nexttile
    y = [i_mean_rest(:); i_mean_task(:)];
    x = [ones(length(i_mean_rest(:)), 1); 2*ones(length(i_mean_task(:)), 1)];
    ylim([0.3 0.75])
    swarmchart(x, y, [], colmat, 'filled')
    xticks([1 2])
    if (i == 6) | (i == 7) | (i == 5), xticklabels(["Rest", "Stimulation"]), else, xticklabels(["", ""]), end
    if (i == 1) | (i == 4) | (i == 7) | (i == 10), ylabel("timescale (\tau)"), end
    stattext = get_stattext_wilcox(zs(i), ps_correct(i), rs(i), 'horizontal');
    ylim([0    0.4])
    sigstar_duodenal({[1 2]}, ps_correct(i), stattext, 0, 0.045);
    % ylim([0    1.5])
    
    title(shufflestring(i))
end
fontsize(gcf, 14, "points")
sgtitle("Recurrent Connections: ")
% saveas(gcf, 'figs\supp_dti_control_for_firing.png')
% saveas(gcf, 'figs\supp_dti_control_for_firing.fig')
exportgraphics(gcf, 'cellreports_review\figs\supp_dti_control_for_firing.jpg')