%% Modeling new figure including both normal and extreme values
clear, clc, close all
cd C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\modeling_3\
addpath C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\
addpath C:\Users\duodenum\Desktop\brain_stuff\misc\slanCM\
normal = load("cellreports_review\dtinetwork.mat");
% more = load("dtinetwork_recurrent.mat");
xtreme = load("cellreports_review\dtinetwork_recurrent_extreme.mat");

nsims=30;
colors = repmat(1:nsims, 360, 1);
colors = [colors(:); colors(:)];
colmat = zeros(length(colors), 3);
all_cols = slanCM('glasbey', 30);

for i = 1:length(colors)
    i_col = colors(i);
    colmat(i, :) = all_cols(i_col, :);
end

acws_rest = permute(xtreme.acwdrs_rest, [1 2 4 3]);
acws_task = permute(xtreme.acwdrs_task, [1 2 4 3]);
%% Calculate mean and sd
var_rest = squeeze(std(acws_rest, [], 2));
var_task = squeeze(std(acws_task, [], 2));

mean_rest = squeeze(mean(acws_rest, 2));
mean_task = squeeze(mean(acws_task, 2));

rt = ((mean_task - mean_rest) ./ mean_rest) * 100;
nshuffles = 11;

%% Do comparisons and corrrection
ps = zeros(1, nshuffles);
rs = zeros(1, nshuffles);
zs = zeros(1, nshuffles);
for i = 1:nshuffles
    i_mean_rest = mean_rest(:, :, i);
    i_mean_task = mean_task(:, :, i);
    [ps(i), ~, stats] = ranksum(i_mean_rest(:), i_mean_task(:));
    zs(i) = stats.zval;
    [~, rs(i)] = ranksum_effect_size(i_mean_rest(:), i_mean_task(:));
end
ps_correct = bonf_holm(ps, 0.05);


%% Make only the first 9

surrogates = [0 0.5 1 1.5 2 2.5 3 3.5 4];
shufflestring = string(surrogates);
close all
figure('Position',[-1.4606e+03 129.8000 784.8000 655.2000]);
plotlocations = [1.5 3.5 6.5 9.5 10.5 11.5 13.5 14.5 15.5];
% subplot(4,4, 'TileSpacing','tight', 'Padding','tight')
textposition = [1.47112462006079,0.290466600026076];
for i = 1:9
    i_mean_rest = mean_rest(:, :, i);
    i_mean_task = mean_task(:, :, i);

    i_rt = rt(:, i, :);
    i_var_rest = var_rest(:, :, i);

    subplot(4, 4, plotlocations(i))
    y = [i_mean_rest(:); i_mean_task(:)];
    x = [ones(length(i_mean_rest(:)), 1); 2*ones(length(i_mean_task(:)), 1)];
    swarmchart(x, y, 8, colmat, 'filled')
    xticks([1 2])
    if (i == 8) | (i == 9) | (i == 7) | (i == 11), xticklabels(["Rest", "Stimulation"]), else, xticklabels(["", ""]), end
    if (i == 1) | (i == 3) | (i == 4) | (i == 7), ylabel("timescale (\tau)"), end
    stattext = get_stattext_wilcox(zs(i), ps_correct(i), rs(i), 'horizontal');
    ylim([0.05    0.3])
    sigstar_duodenal_withpos({[1 2]}, ps_correct(i), stattext, 0, textposition);
    
    title(shufflestring(i))
end
fontsize(gcf, 14, "points")
sgtitle("Recurrent Connections: ")
exportgraphics(gcf, 'cellreports_review\figs\lastfigure.jpg', ...
    'Resolution', 600)
% saveas(gcf, 'figs\lastfigure.fig')

%% Conceptual plot up to 4
rt_2 = mean(rt, [1 2]);
rt_2 = squeeze(rt_2(1, 1, 1:9));
surrogates_2 = surrogates(1:9);
close all
stem(surrogates_2, rt_2, "filled", 'k')
hold on 

plot(0:0.001:4, interp1(surrogates_2 , rt_2, 0:0.001:4, "pchip"), ...
    'Color', "#7E2F8E", 'LineWidth',3)
grid on

xlabel("Strength of Recurrent Connections")
ylabel(["Average % change of timescale (\tau)" "(Stimulation - Rest)"])
fontsize(gcf, 14, "points")
xticks(surrogates_2)
xlim([-0.1 4.1])

exportgraphics(gcf, 'cellreports_review\figs\lastfigure2.jpg', ...
    'Resolution', 600)
% saveas(gcf, 'figs\lastfigure2.fig')