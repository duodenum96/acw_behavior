%% Do SVM for rest / self / other
%% Organize all the data
clear, clc, close all
cd C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\daviddata\rest
rest_acws = load('rest_acws_acwdr');
cd C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\daviddata\task
self_acws = load('self_acws_dr');
other_acws = load('other_acws_dr');
addpath C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\
addpath C:\Users\duodenum\Desktop\brain_stuff\misc\slanCM

rest_all_acws = cat(2, rest_acws.acwdrs_rest{:})';
self_all_acws = cat(2, self_acws.acwdrs_self{:})';
other_all_acws = cat(2, other_acws.acwdrs_other{:})';

nchan = 64;
ntp = 47;
subjcode = repmat((1:21)', ntp, nchan);

nsubj = 21; % a priori
all_acws = [rest_all_acws; self_all_acws; other_all_acws];
subjcode_all = [subjcode; subjcode; subjcode];
labels = [repmat("rest", size(rest_all_acws, 1), 1); 
    repmat("self", size(self_all_acws, 1), 1); 
    repmat("other", size(other_all_acws, 1), 1)];
% subjlabels = [(1:nsubj)'; (1:nsubj)'; (1:nsubj)'];


%% Save data
cd C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\daviddata\
save("acwdata_acwdr.mat")
%% Do Topomaps of ACW values
addpath C:\Users\duodenum\Desktop\brain_stuff\eeglab2023.1
eeglab nogui
addpath('C:\Users\duodenum\Desktop\brain_stuff\eeglab2023.1\functions\sigprocfunc')
self_averaged = mean(self_all_acws);
other_averaged = mean(other_all_acws);
rest_averaged = mean(rest_all_acws);
% Import 1 EEG struct
load('rest\001_rest_reordered.mat')

clims = [min([rest_averaged(:); self_averaged(:); other_averaged(:)]), ...
    max([rest_averaged(:); self_averaged(:); other_averaged(:)])];

close
figure('Position', [180.2000  433.8000  946.4000  305.6000]);
subplot(1, 3, 1)
topoplot(rest_averaged, EEG.chanlocs)
clim(clims)
title('Rest')

subplot(1, 3, 2)
topoplot(self_averaged, EEG.chanlocs)
clim(clims)
title('Self')

subplot(1, 3, 3)
topoplot(other_averaged, EEG.chanlocs)
clim(clims)
title('Other')

colormap summer
cb = colorbar;
ylabel(cb, "timescale (\tau)")
set(cb, 'Position', [ 0.9019    0.2330    0.0156    0.6178])

fontsize(gcf, 14, "points")
% saveas(gcf, 'figs_acwdr\figure1\acw_mean_acwdr.fig')
% saveas(gcf, 'figs_acwdr\figure1\acw_mean_acwdr.jpg')
exportgraphics(gcf, 'figs_acwdr\figure1\acw_mean_acwdr.jpg', 'Resolution',600)

%% ACW change topomaps

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

%% Rest - Other states 

% Correct for multiple comparisons
pvals = zeros(1,3);
[pvals(1),r,U,rsum,h] = ranksum_effect_size(acws_rest_persubj(:), acws_self_persubj(:));
[pvals(2),r,U,rsum,h] = ranksum_effect_size(acws_other_persubj(:), acws_self_persubj(:));
[pvals(3),r,U,rsum,h] = ranksum_effect_size(acws_rest_persubj(:), acws_other_persubj(:));
pvals_correct = bonf_holm(pvals);

allcolors = slanCM('tab20', 21);
col_1subj = repmat((1:21)', 1, nchan);
col_all = [col_1subj(:); col_1subj(:); col_1subj(:)];
colmat_all = zeros(length(col_all), 3);
for i = 1:length(col_all)
    i_subj = col_all(i);
    colmat_all(i, :) = allcolors(i_subj, :);
end

n_acws = length(acws_rest_persubj(:));
data = [acws_rest_persubj(:); acws_self_persubj(:); acws_other_persubj(:)]; 
groups = [ones(n_acws, 1); 2*ones(n_acws, 1); 3*ones(n_acws,1)];
% Save data for mixed effects
memfolder = 'C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\revision\mixed_effects';
channelvar = repmat((1:64), 21, 1); channelvar = repmat(channelvar(:), 3, 1);
tbl = table(data, groups, col_all, channelvar, 'VariableNames', ...
    {'ACW', 'State', 'Subjects', 'Channels'});
writetable(tbl, fullfile(memfolder, 'human_mean_x_time.csv'))

swarmchart(groups, data, [], colmat_all, 'filled')
line([0.75 1.25], [median(data(groups==1)) median(data(groups==1))], ...
    'color', 'k', 'LineWidth', 1.5)
line([1.75 2.25], [median(data(groups==2)) median(data(groups==2))], ...
    'color', 'k', 'LineWidth', 1.5)
line([2.75 3.25], [median(data(groups==3)) median(data(groups==3))], ...
    'color', 'k', 'LineWidth', 1.5)
set(gca, 'YGrid', 'on')

% 1-2
[~, ~, stats] = ranksum(acws_rest_persubj(:), acws_self_persubj(:));
[~,r,U,rsum,h] = ranksum_effect_size(acws_rest_persubj(:), acws_self_persubj(:));
pval = pvals_correct(1);
stattext = get_stattext_wilcox(stats.zval, pval, r, 'horizontal');
sigstar_duodenal({[1 2]}, pval,stattext, 0)
% 2-3
[~, ~, stats] = ranksum(acws_other_persubj(:), acws_self_persubj(:));
[~,r,U,rsum,h] = ranksum_effect_size(acws_other_persubj(:), acws_self_persubj(:));
pval = pvals_correct(2);
stattext = get_stattext_wilcox(stats.zval, pval, r, 'horizontal');
sigstar_duodenal({[2 3]}, pval,stattext, 0)
% 1-3
[~, ~, stats] = ranksum(acws_rest_persubj(:), acws_other_persubj(:));
[~,r,U,rsum,h] = ranksum_effect_size(acws_rest_persubj(:), acws_other_persubj(:));
pval = pvals_correct(3);
stattext = get_stattext_wilcox(stats.zval, pval, r, 'horizontal');
sigstar_duodenal({[1 3]}, pval,stattext, 0)

% kruskalwallis(data, groups)

xticks([1 2 3])
xticklabels(["Rest", "Self", "Other"])
ylabel("Mean of timescale (\tau) across windows")

fontsize(gcf, 14, "points")
% saveas(gcf, 'figs_acwdr\figure1\acwdiff_stat_acwdr.fig')
% saveas(gcf, 'figs_acwdr\figure1\acwdiff_stat_acwdr.jpg')
exportgraphics(gcf, 'figs_acwdr\figure1\acwdiff_stat_acwdr.jpg', 'Resolution', 600)
%% Same as above but for average across scalp
allcolors = slanCM('glasbey', 21);
col_1subj = repmat((1:21)', 1, ntp);
col_all = [col_1subj(:); col_1subj(:); col_1subj(:)];
colmat_all = zeros(length(col_all), 3);
for i = 1:length(col_all)
    i_subj = col_all(i);
    colmat_all(i, :) = allcolors(i_subj, :);
end

% Correct for multiple comparisons
pvals = zeros(1,3);
[pvals(1),r,U,rsum,h] = ranksum_effect_size(acws_rest_x_scalp(:), acws_self_x_scalp(:));
[pvals(2),r,U,rsum,h] = ranksum_effect_size(acws_other_x_scalp(:), acws_self_x_scalp(:));
[pvals(3),r,U,rsum,h] = ranksum_effect_size(acws_rest_x_scalp(:), acws_other_x_scalp(:));
pvals_correct = bonf_holm(pvals);

n_acws = length(acws_rest_x_scalp(:));
data = [acws_rest_x_scalp(:); acws_self_x_scalp(:); acws_other_x_scalp(:)]; 
groups = [ones(n_acws, 1); 2*ones(n_acws, 1); 3*ones(n_acws,1)];
% Save data for mixed effects
memfolder = 'C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\revision\mixed_effects';
channelvar = repmat((1:64), 21, 1); channelvar = repmat(channelvar(:), 3, 1);
tbl = table(data, groups, col_all, 'VariableNames', ...
    {'ACW', 'State', 'Subjects'});
writetable(tbl, fullfile(memfolder, 'human_mean_x_space.csv'))

swarmchart(groups, data, [], colmat_all, 'filled')

line([0.75 1.25], [median(data(groups==1)) median(data(groups==1))], ...
    'color', 'k', 'LineWidth', 1.5)
line([1.75 2.25], [median(data(groups==2)) median(data(groups==2))], ...
    'color', 'k', 'LineWidth', 1.5)
line([2.75 3.25], [median(data(groups==3)) median(data(groups==3))], ...
    'color', 'k', 'LineWidth', 1.5)
set(gca, 'YGrid', 'on')

% 1-2
[~, ~, stats] = ranksum(acws_rest_x_scalp(:), acws_self_x_scalp(:));
[~,r,U,rsum,h] = ranksum_effect_size(acws_rest_x_scalp(:), acws_self_x_scalp(:));
pval = pvals_correct(1);
stattext = get_stattext_wilcox(stats.zval, pval, r, 'horizontal');
sigstar_duodenal({[1 2]}, pval,stattext, 0)
% 2-3
[~, ~, stats] = ranksum(acws_self_x_scalp(:), acws_other_x_scalp(:));
[~,r,U,rsum,h] = ranksum_effect_size(acws_self_x_scalp(:), acws_other_x_scalp(:));
pval = pvals_correct(2);
stattext = get_stattext_wilcox(stats.zval, pval, r, 'horizontal');
sigstar_duodenal({[2 3]}, pval,stattext, 0)
% 1-3
[~, ~, stats] = ranksum(acws_rest_x_scalp(:), acws_other_x_scalp(:));
[~,r,U,rsum,h] = ranksum_effect_size(acws_rest_x_scalp(:), acws_other_x_scalp(:));
pval = pvals_correct(3);
stattext = get_stattext_wilcox(stats.zval, pval, r, 'horizontal');
sigstar_duodenal({[1 3]}, pval,stattext, 0)

xticks([1 2 3])
xticklabels(["Rest", "Self", "Other"])
ylabel("Mean of timescale (\tau) across channels")

fontsize(gcf, 14, "points")
% saveas(gcf, 'figs_acwdr\figure1\acwdiff_x_scalp_stat_acwdr.fig')
% saveas(gcf, 'figs_acwdr\figure1\acwdiff_x_scalp_stat_acwdr.jpg')
exportgraphics(gcf, 'figs_acwdr\figure1\acwdiff_x_scalp_stat_acwdr.jpg', ...
    'Resolution',600)

kruskalwallis(data, groups)
% chi^2 = 77.56, p<0.001

%% Rest SD - rest-task change correlation

allcolors = slanCM('tab20', nsubj);
col_1subj = repmat((1:21)', 1, nchan);
col_all = col_1subj(:);
% col_all = [col_1subj(:); col_1subj(:); col_1subj(:)];
colmat_all = zeros(length(col_all()), 3);
for i = 1:length(col_all)
    i_subj = col_all(i);
    colmat_all(i, :) = allcolors(i_subj, :);
end

% correct p values
pvals = zeros(2, 1);

[~, pvals(1)] = corr(acws_sd_rest_persubj(:), restself_mean(:),'type','Spearman');
[~, pvals(2)] = corr(acws_sd_rest_persubj(:), restother_mean(:),'type','Spearman');
pvals_correct = bonf_holm(pvals);


% Variables: acws_sd_rest_persubj, restself_mean, restother_mean
close all
figure;
subplot(1,2, 1)
x = acws_sd_rest_persubj(:);
y = restself_mean(:);
scatter(x, y, [], colmat_all, 'filled')
hold on

[no_outlier_data, out_idx] = ...
        rmoutliers([acws_sd_rest_persubj(:), restself_mean(:)]);
scatter(x(out_idx), y(out_idx), 'ko')
xlabel('Rest STD of timescale (\tau)')
ylabel({'% Change of timescale (\tau)', '(Task - Rest)'})
axis square
title('Self - Rest')
[rho, ~] = corr(acws_sd_rest_persubj(:), restself_mean(:),'type','Spearman');
pval = pvals_correct(1);
assert(pval < 0.001) % if not, fix the sprintf
xticks(0:10:50)
stattext = sprintf('$ \\rho = %.2f, p < 0.001$', rho);
text(min(xlim) + 0.05, min(ylim)+0.05, stattext, ...
            'Horiz','left', 'Vert','bottom', 'Interpreter', 'LaTeX')
grid on
hold on
rl = refline(polyfit(x, y, 1));
set(rl, 'Color', 'k')
xticks(0:0.02:0.1)

subplot(1,2, 2)
x = acws_sd_rest_persubj(:);
y = restother_mean(:);
[no_outlier_data, out_idx] = ...
        rmoutliers([acws_sd_rest_persubj(:), restother_mean(:)]);
scatter(x, y, [], colmat_all, 'filled')
hold on
scatter(x(out_idx), y(out_idx), 'ko')
xlabel('Rest STD of timescale (\tau)')
axis square
title('Other - Rest')
[rho, ~] = corr(acws_sd_rest_persubj(:), restother_mean(:),'type','Spearman');
pval = pvals_correct(2);
assert(pval < 0.001) % if not, fix the sprintf
xticks(0:10:50)
stattext = sprintf('$ \\rho = %.2f, p < 0.001$', rho);
text(min(xlim) + 0.05, min(ylim)+0.05, stattext, ...
            'Horiz','left', 'Vert','bottom', 'Interpreter', 'LaTeX')
grid on
hold on
rl = refline(polyfit(x, y, 1));
set(rl, 'Color', 'k')
xticks(0:0.02:0.1)


fontsize(gcf, 14, "points")
% saveas(gcf, 'figs_acwdr\figure1\restacw_resttaskchange_acwdr.fig')
% saveas(gcf, 'figs_acwdr\figure1\restacw_resttaskchange_acwdr.jpg')
exportgraphics(gcf, 'figs_acwdr\figure1\restacw_resttaskchange_acwdr.jpg', ...
    'Resolution',600)
%% Without outliers

% correct p values
pvals = zeros(1, 2);
[no_outlier_data, out_idx] = ...
        rmoutliers([acws_sd_rest_persubj(:), restself_mean(:)]);
[~, pvals(1)] = corr(no_outlier_data(:, 1), no_outlier_data(:,2),'type','Spearman');
[no_outlier_data, out_idx] = ...
        rmoutliers([acws_sd_rest_persubj(:), restother_mean(:)]);
[~, pvals(2)] = corr(no_outlier_data(:, 1), no_outlier_data(:,2),'type','Spearman');
pvals_correct = bonf_holm(pvals);


close all
figure;
subplot(1,2, 1)
[no_outlier_data, out_idx] = ...
        rmoutliers([acws_sd_rest_persubj(:), restself_mean(:)]);
[rho1, ~] = corr(no_outlier_data(:, 1), no_outlier_data(:,2),'type','Spearman');
pval1 = pvals_correct(1);
scatter(no_outlier_data(:, 1), no_outlier_data(:, 2), [], colmat_all(~out_idx, :), 'filled')
xlabel('Rest STD of timescale (\tau)')
ylabel({'% Change of timescale (\tau)', '(Task - Rest)'})
axis square
title('Self - Rest')
if pval1 < 0.001
    stattext = sprintf('$ \\rho = %.2f, p < 0.001$', rho1);
else
    stattext = sprintf('$ \\rho = %.2f, p = %.3f $', rho1, pval1);
end
text(min(xlim) + 0.05, min(ylim)+0.05, stattext, ...
            'Horiz','left', 'Vert','bottom', 'Interpreter', 'LaTeX')
grid on
hold on
rl = refline(polyfit(no_outlier_data(:, 1), no_outlier_data(:, 2), 1));
set(rl, 'Color', 'k')


subplot(1,2, 2)
[no_outlier_data, out_idx] = ...
        rmoutliers([acws_sd_rest_persubj(:), restother_mean(:)]);
[rho2, ~] = corr(no_outlier_data(:, 1), no_outlier_data(:,2),'type','Spearman');
pval2 = pvals_correct(2);
scatter(no_outlier_data(:, 1), no_outlier_data(:, 2), [], colmat_all(~out_idx, :), 'filled')
xlabel('Rest STD of timescale (\tau)')
axis square
title('Other - Rest')
if pval2 < 0.001
    stattext = sprintf('$ \\rho = %.2f, p < 0.001$', rho2);
else
    stattext = sprintf('$ \\rho = %.2f, p = %.3f $', rho2, pval2);
end
text(min(xlim) + 0.05, min(ylim)+0.05, stattext, ...
            'Horiz','left', 'Vert','bottom', 'Interpreter', 'LaTeX')
grid on
hold on
rl = refline(polyfit(no_outlier_data(:, 1), no_outlier_data(:, 2), 1));
set(rl, 'Color', 'k')

fontsize(gcf, 14, "points")

% saveas(gcf, 'figs_acwdr\figure1\restacw_resttaskchange_acwdr_no_outliers.fig')
% saveas(gcf, 'figs_acwdr\figure1\restacw_resttaskchange_acwdr_no_outliers.jpg')
exportgraphics(gcf, 'figs_acwdr\figure1\restacw_resttaskchange_acwdr_no_outliers.jpg', ...
    'Resolution',600)



%% Topomap of mean acw change and rest acw sd

restself_med = median(restself_mean);
restother_med = median(restother_mean);

rest_sd_med = median(acws_sd_rest_persubj);

clims = [min([restself_med(:); restother_med(:)]), ...
    max([restself_med(:); restother_med(:)])];


close
figure('Position', [180.2000  433.8000  946.4000  305.6000]);
ax1 = subplot(1, 3, 1);
topoplot(rest_sd_med, EEG.chanlocs)
colorbar
clim([min(rest_sd_med) max(rest_sd_med)])
title('Rest STD of timescale (\tau)')
colormap(ax1, 'summer')

ax2 = subplot(1, 3, 2);
topoplot(restself_med, EEG.chanlocs)
clim([-max(abs(clims)), max(abs(clims))])
title('Self - Rest')
colormap(ax2, slanCM('vik'))

ax3 = subplot(1, 3, 3);
topoplot(restother_med, EEG.chanlocs)
clim([-max(abs(clims)), max(abs(clims))])
title('Other - Rest')
colormap(ax3, slanCM('vik'))
colormap(ax2, slanCM('vik'))
colormap(ax1, 'summer')

cb = colorbar;
ylabel(cb, ["% Change of timescale (\tau)",  "(Task - Rest)"])
set(cb, 'Position', [ 0.9019    0.2330    0.0156    0.6178])
set(cb, 'Limit', clims)

fontsize(gcf, 14, "points")
% saveas(gcf, 'figs_acwdr\figure1\acw_std_prcchange_median.fig')
% saveas(gcf, 'figs_acwdr\figure1\acw_std_prcchange_median.jpg')
exportgraphics(gcf, 'figs_acwdr\figure1\acw_std_prcchange_median.jpg', ...
    'Resolution', 600) 