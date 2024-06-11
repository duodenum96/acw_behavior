%% Analyze data from DTI network
clear, clc, close all
cd C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\modeling_3
addpath C:\Users\duodenum\Desktop\brain_stuff\misc
addpath C:\Users\duodenum\Desktop\brain_stuff\misc\slanCM
addpath 'C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw'
nsims = 30;
colors = repmat(1:nsims, 360, 1);
colors = [colors(:); colors(:)];
colmat = zeros(length(colors), 3);
all_cols = slanCM('glasbey', nsims);

for i = 1:length(colors)
    i_col = colors(i);
    colmat(i, :) = all_cols(i_col, :);
end


load cellreports_review\dtinetwork.mat

var_rest = squeeze(std(acwdrs_rest, [], 2));
var_task = squeeze(std(acwdrs_task, [], 2));

mean_rest = squeeze(mean(acwdrs_rest, 2));
mean_task = squeeze(mean(acwdrs_task, 2));

rtchange = ((mean_task - mean_rest) ./ mean_rest) .* 100;
%% Plot mean difference
close all
y = [mean_rest(:); mean_task(:)];
x = [ones(length(mean_rest(:)), 1); 2*ones(length(mean_task(:)), 1)];

swarmchart(x, y, [], colmat, 'filled')
xticks([1 2])
xticklabels(["Rest", "Stimulation"])
ylabel("timescale (\tau)")
% ylim([0.06 0.12])

[p, ~, stats] = ranksum(mean_rest(:), mean_task(:));
[~, r] = ranksum_effect_size(mean_rest(:), mean_task(:));
stattext = get_stattext_wilcox(stats.zval, p, r, 'horizontal');
sigstar_duodenal({[1 2]}, p, stattext, 0);
fontsize(gcf, 14, "points")
% ylim([0.07 0.12])
title("Recurrent Connections: 1")
saveas(gcf, 'cellreports_review\figs\dti_healthy.jpg')
% saveas(gcf, 'cellreports_review\figs\dti_healthy.fig')
% ylim: 0.3 - 0.75
exportgraphics(gcf, 'cellreports_review\figs\dti_healthy.jpg', 'Resolution',600)
%% Scatter mean
close all
[rho, pval] = corr(var_rest(:), rtchange(:));
scatter(var_rest(:), rtchange(:), [], colmat(1:10800, :), 'filled')
% xlim([0 0.8])
grid on
line = lsline;
set(line, 'color', 'r')
stattext = sprintf("$ \\rho = %.3f, p < 0.001 $", rho);
text(min(xlim) + 0.001, min(ylim)+0.01, stattext, ...
            'Horiz','left', 'Vert','bottom', 'Interpreter', 'LaTeX')
fontsize(gcf, 14, "points")
xlabel("STD of timescale (\tau) across windows during rest")
ylabel(["% change of timescale (\tau)" "(Stimulation - Rest)"])
% saveas(gcf, 'figs\dti_healthy_corr.jpg')
% saveas(gcf, 'figs\dti_healthy_corr.fig')
xlim([0, 0.13])
ylim([-25 70])
exportgraphics(gcf, 'cellreports_review\figs\dti_healthy_corr.jpg', 'Resolution', 600)
% xlim: 0-0.8
% ylim: -40 100