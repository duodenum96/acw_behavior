%% Pipeline for all the figures
%% ML figure
clear, clc, close all
cd C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw
load("acw0_ca.mat")

labels = ["onset", "offset", "initial rest", "sustained rest", "locomotion"];
labels_realdeal = strings(size(labels_onehot, 1), 1);
for i = 1:length(labels)
    indx = labels_onehot == i;
    labels_realdeal(indx) = labels(i);
end

%% Supplementary material: topomaps

mice = {'cm124', 'cm125', 'cm126', 'cm127', 'cm128'};
nmice = length(mice);
topo = cell(1, nmice);

for i = 1:nmice
    files = dir("C:\Users\duodenum\Desktop\brain_stuff\mice_wfoi\WFOM_extracted_signals\WFOM_extracted_signals\" + mice{i} + "*");
    data = load("C:\Users\duodenum\Desktop\brain_stuff\mice_wfoi\WFOM_extracted_signals\WFOM_extracted_signals\" + files(1).name);
    topo{i} = data.info.rois;
end

roidims = size(topo{1});
nroi = 92;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% average topographies for each state %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nlabels = length(labels);
% acwtopos = zeros([roidims, nmice, nlabels]);
acwtopos = cell(1, nmice);
for i = 1:nmice
    mouse = mice{i};
    mouse_indx = string(mousenames) == string(mouse);
    acwtopos{i} = cell(1, nlabels);
    for j = 1:nlabels
        acwtopos{i}{j} = zeros(roidims);
        label_indx = labels_onehot == j;
        acwtopo = mean(acwmat(mouse_indx & label_indx, :));
        for k = 1:nroi
            roi_indx = topo{i} == k;
            acwtopos{i}{j}(roi_indx) = acwtopo(k);
        end
    end
end
% Set nonrois to nan
for i = 1:nmice
    for j = 1:nlabels
        indx = acwtopos{i}{j} == 0;
        acwtopos{i}{j}(indx) = nan;
    end
end
%% Make nice plots
js = [1 5 2 3 4]; % onset, run, offset, initrest, susrest

close all
for i = 1:nmice
    figure(Position= [-1588 378 1320 420]);
    tiledlayout(1, nlabels, "TileSpacing","tight", "Padding", "tight")
    for j  = js
        nexttile
        colormap summer
        h = imagesc(acwtopos{i}{j});
        set(h, 'AlphaData', 1 - isnan(acwtopos{i}{j})) % get rid of background
        clim([min([acwtopos{i}{:}], [],  "all"), max([acwtopos{i}{:}], [],  "all")])
        if j == js(end), cb = colorbar; ylabel(cb, 'ACW-0 (s)'), end
        axis square
        axis off
        title(labels(j))
    end
    sgtitle(mice{i})
    fontsize(gcf, 13, "points")
    % saveas(gcf, ['figs\figs_clean_all\topomap_', mice{i}, '.png'])
    % saveas(gcf, ['figs\figs_clean_all\topomap_', mice{i}, '.fig'])
    close
end
% save("topo_acw0.mat") % save intermediate results
%% Calculate variabilities
% Prepare the data
addpath('C:\Users\duodenum\Desktop\brain_stuff\misc\hex_and_rgb')
load('yaxis_info.mat', 'order')
sessions = string(sessions);
labels = ["onset", "offset", "initial rest", "sustained rest", "running"];

labels_names = strings(size(labels_onehot, 1), 1);
for i = 1:length(labels)
    indx = labels_onehot == i;
    labels_names(indx) = labels(i);
end

dict = table(labels_names, acwmat, string(mousenames), sessions, ...
    'VariableNames', {'labels', 'acw', 'mice', 'sessions'});
%% Calculate rest variability
susrest = dict.labels == "sustained rest";
run = dict.labels == "running";
sesh = dict.sessions;
% Rest variability x trials for each ROI
rest = dict(susrest, :);
restacws = rest.acw;
restmice = rest.mice;
restses = sesh(susrest);

rundata = dict(run, :);
runacws = rundata.acw;
runmice = rundata.mice;
runses = sesh(run);


% For every mouse, calculate variability of every ROI across trials
acw_var_x_trials = zeros(nroi, nmice);
acw_mean_x_trials = zeros(nroi, nmice);
acw_mean_x_trials_running = zeros(nroi, nmice);
for i = 1:nmice
    mouse = string(mice{i});
    mouseindx = restmice == mouse;
    acws = restacws(mouseindx, :);
    acw_var_x_trials(:, i) = std(acws);
    acw_mean_x_trials(:, i) = mean(acws);

    mouseindx = runmice == mouse;
    acws = runacws(mouseindx, :);
    acw_mean_x_trials_running(:, i) = mean(acws);
end
% old name: mice_onehot_trial
mice_onehot_x_trial = [ones(nroi, 1); 2*ones(nroi, 1); 3*ones(nroi, 1); 4*ones(nroi, 1); 5*ones(nroi, 1)];

% Same for across ROIs
acw_var_x_rois = cell(nmice, 1);
acw_mean_x_rois = cell(nmice, 1);
acw_mean_x_rois_running = cell(nmice, 1);
ntrials = zeros(nmice, 1);
ntrials_running = zeros(nmice, 1);
for i = 1:nmice
    mouse = string(mice{i});
    mouseindx = restmice == mouse;
    acws = restacws(mouseindx, :);
    acw_var_x_rois{i} = std(acws, [], 2);
    acw_mean_x_rois{i} = mean(acws, 2);
    ntrials(i) = length(acw_var_x_rois{i});

    mouseindx = runmice == mouse;
    acws = runacws(mouseindx, :);
    acw_mean_x_rois_running{i} = mean(acws, 2);
end
% old name: mice_onehot
mice_onehot_x_roi = [ones(ntrials(1), 1); 2*ones(ntrials(2), 1); ...
    3*ones(ntrials(3), 1); 4*ones(ntrials(4), 1); 5*ones(ntrials(5), 1)];

mice_onehot_x_roi_run = [ones(ntrials_running(1), 1); 2*ones(ntrials_running(2), 1); ...
    3*ones(ntrials_running(3), 1); 4*ones(ntrials_running(4), 1); 5*ones(ntrials_running(5), 1)];

 % In the end, we have acw_var_x_rois, acw_mean_x_rois, acw_var_x_trials,
 % acw_mean_x_trial and the *_running variants of these

%% Get colors
% This is the mice indices for x_trials
colors = [102,194,165
          252,141,98
          141,160,203
          231,138,195
          166,216,84] ./ 255;
ncol = size(mice_onehot_x_trial, 1);
colmat = zeros(ncol, 3);
for i = 1:nmice
    colmat(mice_onehot_x_trial == i, :) = repmat(colors(i, :), nroi, 1);
end
%% save
calculations = [];
calculations.acw_var_x_trials = acw_var_x_trials;
calculations.acw_var_x_rois = acw_var_x_rois;
calculations.acw_var_x_rois_col = cat(1, acw_var_x_rois{:});
calculations.acw_var_x_rois_mice = mice_onehot_x_roi;
calculations.acw_var_x_trials_mice = mice_onehot_x_trial;
calculations.acw_var_x_trials_col = acw_var_x_trials(:);
calculations.acw_mean_x_rois_col = cat(1, acw_mean_x_rois{:});
calculations.acw_mean_x_rois = acw_mean_x_rois;

% save('restvariability.mat', "calculations")

calculations = [];
calculations.acw_mean_x_rois_running = acw_mean_x_rois_running;
calculations.acw_mean_x_trials_running = acw_mean_x_trials_running;
calculations.acw_mean_x_rois_running_col = cat(1, acw_mean_x_rois_running{:});
calculations.mice_onehot_x_roi_run = mice_onehot_x_roi_run;

% save('runvariability.mat', "calculations")

%% Calculate prc change of mean 
acw_mean_prc_x_trl = ((acw_mean_x_trials_running - acw_mean_x_trials) ./ ...
    acw_mean_x_trials) .* 100;
%% Correlate rest variability to rest-task change (MEANS)
addpath('C:\Users\duodenum\Desktop\brain_stuff\misc\slanCM')
h = scatterhist(acw_var_x_trials(:), acw_mean_prc_x_trl(:), ...
    'group', mice_onehot_x_trial, 'kernel', 'on', ...
    'Direction', 'out');
hold on
rl = refline(polyfit(acw_var_x_trials(:), acw_mean_prc_x_trl(:), 1));
set(rl, 'Color', 'k')
% mdl = fitlm(acw_var_x_trials(:), acw_mean_prc_x_trl(:));
% hold on
% plot(mdl)

xlabel('STD of ACW-0 across windows during rest')
ylabel(["% ACW-0 change" "from rest to locomotion"])

scs = get(h(1), 'Children'); scs = scs(2:end);

for i = 1:5
    set(scs(i), "MarkerFaceColor", colors(i, :));
    set(scs(i), "MarkerEdgeColor", "none");
    set(h(2).Children(i+1), "Color", colors(i, :))
    set(h(3).Children(i+1), "Color", colors(i, :))
end
set(h(1).Legend, 'String', mice, 'box', 'off')
grid on

[rho, pval] = corr(acw_var_x_trials(:), acw_mean_prc_x_trl(:), 'Type', 'Spearman');
text(min(xlim) + 0.05, min(ylim)+0.05, sprintf('$ \\rho = %.2f, p < 0.001$', rho), ...
            'Horiz','left', 'Vert','bottom', 'Interpreter', 'LaTeX')
fontsize(gcf, 13, "points")
saveas(gcf, 'figs\figs_clean_all\restrun_var_prcchange_corr.fig')
saveas(gcf, 'figs\figs_clean_all\restrun_var_prcchange_corr.png')

%% Same for individual mice (SD - Mean Change)
close all
figure('Position', [189 398 1314 360])
tiledlayout(1, nmice)
for i = 1:nmice
    nexttile
    scatter(acw_var_x_trials(:, i), acw_mean_prc_x_trl(:, i), ...
        [], colors(i, :), 'filled')
    rl = refline(polyfit(acw_var_x_trials(:, i), acw_mean_prc_x_trl(:, i), 1));
    set(rl, 'Color', 'k')

    axis square
    [rho, pval] = corr(acw_var_x_trials(:, i), acw_mean_prc_x_trl(:, i));
    assert(pval < 0.001, "p value not under 0.001, fix line 250")
    text(min(xlim) + 0.05, min(ylim), ...
        sprintf('$ \\rho = %.2f, p < 0.001$', rho), ...
                'Horiz','left', 'Vert','bottom', 'Interpreter', 'LaTeX') % bonf. corrected
    if i == 3, xlabel('STD of ACW-0 across windows during rest'), end
    if i == 1, ylabel({'% ACW-0 change', 'from rest to locomotion'}), end
    grid on
    title(mice(i))
end
fontsize(gcf, 13, "points")
saveas(gcf, 'figs\figs_clean_all\restrun_var_prcchange_corr_individualmice.fig')
saveas(gcf, 'figs\figs_clean_all\restrun_var_prcchange_corr_individualmice.png')

%% Get MORE colors (because why not)
colors_sd = slanCM('Set3', 5);

ncol = size(mice_onehot_x_trial, 1);
colmat_sd = zeros(ncol, 3);
for i = 1:nmice
    colmat_sd(mice_onehot_x_trial == i, :) = repmat(colors_sd(i, :), nroi, 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rest - Running change plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('acw0_ca.mat', 'sessions')
runtrl = labels_realdeal == "locomotion";
resttrl = labels_realdeal == "sustained rest";

mousenames = string(mousenames);
acw_run = acwmat(runtrl, :);
mousenames_run = mousenames(runtrl);
acw_rest = acwmat(resttrl, :);
mousenames_rest = mousenames(resttrl);

mousenames_run_code = zeros(length(mousenames_run), nroi);
for i = 1:nmice, mousenames_run_code(mousenames_run == string(mice(i)), :) = i; end
mousenames_rest_code = zeros(length(mousenames_rest), nroi);
for i = 1:nmice, mousenames_rest_code(mousenames_rest == string(mice(i)), :) = i; end

%% Average across trials

sessions = string(sessions);
restsessions = sessions(resttrl);
uniqtrials = unique(restsessions);
n_trials = length(uniqtrials);
rest_averages = zeros(n_trials, nroi);
mousenames_averages_rest =  zeros(n_trials, 1);

for i =1:n_trials
    indx = restsessions == uniqtrials(i);
    tmp = mousenames_rest_code(indx, 1);
    mousenames_averages_rest(i) = tmp(1);
    rest_averages(i,:) = mean(acw_rest(indx, :), 1); % For each ROI, average across trials
end

runsessions = sessions(runtrl);
uniqtrials = unique(runsessions);
n_trials = length(uniqtrials);
run_averages = zeros(n_trials, nroi);
mousenames_averages_run =  zeros(n_trials, 1);
for i =1:n_trials
    indx = runsessions == uniqtrials(i);
    tmp = mousenames_run_code(indx, 1);
    mousenames_averages_run(i) = tmp(1);
    run_averages(i, :) = mean(acw_run(indx, :), 1); % For each ROI, average across trials
end
%% Average across trials Rest Run diff
nrois = 92;
close all
restacw_ave = zeros(nrois, nmice);
runacw_ave = zeros(nrois, nmice);
for i = 1:nmice
    i_restmice = mousenames_rest_code(:, 1) == i;
    i_runmice = mousenames_run_code(:, 1) == i;
    restacw_ave(:, i) = mean(acw_rest(i_restmice, :), 1);
    runacw_ave(:, i) = mean(acw_run(i_runmice, :), 1);
end
mice_code = repmat(1:5, nrois, 1);

y = [restacw_ave(:); runacw_ave(:)];
x = [ones(length(restacw_ave(:)), 1); 2 * ones(length(restacw_ave(:)), 1)];
mice_code2 = [mice_code(:); mice_code(:)];

[~, ~, stats] = ranksum(restacw_ave(:), runacw_ave(:));
[pval,r,U,rsum,h] = ranksum_effect_size(restacw_ave(:), runacw_ave(:));
stattext = get_stattext_wilcox(stats.zval, pval, r, 'horizontal');
% Get colmat
colmat = zeros(length(x), 3);
micecol = mice_code2(:);
for i = 1:length(x)
    i_mice = micecol(i);
    colmat(i, :) = colors(i_mice, :);
end
figure(Position=[488 266 739 492])
swarmchart(x, y, [], colmat, 'filled')
xticks([1 2])
xticklabels(["Sustained Rest", "Locomotion"])
sigstar_duodenal({[1 2]}, pval,stattext, 0)
ylabel("Mean of ACW-0 across windows (s)")
hold on
line([0.75 1.25], [median(y(x==1)) median(y(x==1))], 'color', 'k', 'LineWidth', 1.5)
line([1.75 2.25], [median(y(x==2)) median(y(x==2))], 'color', 'k', 'LineWidth', 1.5)
set(gca, 'YGrid', 'on')
legend_sw(colors)
% title(mice(i))
% xtickangle(15)

fontsize(gcf, 14, "points")
saveas(gcf, "figs\figs_clean_all\mean_x_trl_allmice.png")
saveas(gcf, "figs\figs_clean_all\mean_x_trl_allmice.fig")
%% For each mouse, average across trials, rest run difference
%% Do correction for multiple comparisons
pvals = zeros(1, nmice);
for i = 1:nmice
    i_restmice = mousenames_rest_code(:, 1) == i;
    i_runmice = mousenames_run_code(:, 1) == i;
    i_restacw = acw_rest(i_restmice, :);
    i_runacw = acw_run(i_runmice, :);
    [pvals(i),~,~,~,~] = ranksum_effect_size(mean(i_restacw, 1)', mean(i_runacw, 1)');
end
pvals_correct = bonf_holm(pvals);

%% Plot
close all
figure('Position',[81.8000 398.6000 1352 359.4000])
tiledlayout(1, nmice)
for i = 1:nmice
    nexttile
    pval = pvals_correct(i);
    i_restmice = mousenames_rest_code(:, 1) == i;
    i_runmice = mousenames_run_code(:, 1) == i;
    i_restacw = acw_rest(i_restmice, :);
    i_runacw = acw_run(i_runmice, :);

    y = [mean(i_restacw, 1)'; mean(i_runacw, 1)'];
    x = [ones(nroi, 1); 2 * ones(nroi, 1)];

    [~, ~, stats] = ranksum(mean(i_restacw, 1)', mean(i_runacw, 1)');
    [~,r,U,rsum,h] = ranksum_effect_size(mean(i_restacw, 1)', mean(i_runacw, 1)');
    stattext = get_stattext_wilcox(stats.zval, pval, r, 'vertical');

    swarmchart(x, y, [], colors(i, :), 'filled')
    xticks([1 2])
    xticklabels(["Sustained Rest", "Locomotion"])
    sigstar_duodenal({[1 2]}, pval,stattext, 0)
    if i == 1, ylabel("Mean of ACW-0 across windows (s)"), end
    % title(mice(i))
    xtickangle(15)
    line([0.75 1.25], [median(y(x==1)) median(y(x==1))], 'color', 'k', 'LineWidth', 1.5)
    line([1.75 2.25], [median(y(x==2)) median(y(x==2))], 'color', 'k', 'LineWidth', 1.5)
    set(gca, 'YGrid', 'on')
end
fontsize(gcf, 14, "points")
saveas(gcf, "figs\figs_clean_all\rest_run_comparison_meantrials_indmice.png")
saveas(gcf, "figs\figs_clean_all\rest_run_comparison_meantrials_indmice.fig")

%% Mean across ROIs
% Get colors

y = [restmean_x_roi; runmean_x_roi];
x = [ones(length(restmean_x_roi), 1); 2 * ones(length(runmean_x_roi), 1)];
mice_i = [restvar_x_roi_mice; runvar_x_roi_mice];
% colmat = colorgenerator([], mice_i);
colors_4 = slanCM('vivid', 5);

ncol = size(x, 1);
colmat_4 = zeros(ncol, 3);
for i = 1:length(x)
    i_mouse = mice_i(i);
    colmat_4(i, :) = colors_4(i_mouse, :);
end


figure(Position=[488 266 739 492])
swarmchart(x, y, [], colmat_4, 'filled')
xticks([1 2])
xticklabels(["Sustained Rest", "Locomotion"])
ylabel("Mean of ACW-0 across ROIs (s)")
[p, ~, stats] = ranksum(restmean_x_roi, runmean_x_roi);
[~, r] = ranksum_effect_size(restmean_x_roi, runmean_x_roi);
stattext = get_stattext_wilcox(stats.zval, p, r, 'horizontal');
sigstar_duodenal({[1 2]}, p, stattext, 0)


hold on
line([0.75 1.25], [median(y(x==1)) median(y(x==1))], 'color', 'k', 'LineWidth', 1.5)
line([1.75 2.25], [median(y(x==2)) median(y(x==2))], 'color', 'k', 'LineWidth', 1.5)
set(gca, 'YGrid', 'on')

legend_sw(colors_4);
fontsize(gcf, 14, "points")
saveas(gcf, "figs\figs_clean_all\rest_run_mean_roi.png")
saveas(gcf, "figs\figs_clean_all\rest_run_mean_roi.fig")
%% Do for individual mice (Mean x ROIs)
ps = zeros(1, nmice);
rs = zeros(1, nmice);
zs = zeros(1, nmice);
for i = 1:nmice
    mice_ii = mice_i == i;
    y_i = y(mice_ii);
    x_i = x(mice_ii); % x_i == 1: rest, x_i == 2: locomotion
    [ps(i), ~, stats] = ranksum(y_i(x_i==1), y_i(x_i==2));
    zs(i) = stats.zval;
    [~, rs(i)] = ranksum_effect_size(y_i(x_i==1), y_i(x_i==2));
end

ps_correct = bonf_holm(ps, 0.05);

close all
figure('Position',[81.8000 398.6000 1352 359.4000])
tiledlayout(1, nmice)
for i = 1:nmice
    nexttile
    mice_ii = mice_i == i;
    y_i = y(mice_ii);
    x_i = x(mice_ii);

    swarmchart(x_i, y_i, [], colors_4(i, :), 'filled')
    xticks([1 2])
    xticklabels(["Sustained Rest", "Locomotion"])
    if i == 1, ylabel("Mean of ACW-0 across ROIs (s)"), end
    stattext = get_stattext_wilcox(zs(i), ps(i), rs(i), 'vertical');

    hold on
    line([0.75 1.25], [median(y_i(x_i==1)) median(y_i(x_i==1))], 'color', 'k', 'LineWidth', 1.5)
    line([1.75 2.25], [median(y_i(x_i==2)) median(y_i(x_i==2))], 'color', 'k', 'LineWidth', 1.5)
    set(gca, 'YGrid', 'on')

    sigstar_duodenal({[1 2]}, ps_correct(i), stattext, 0)
    % if i == 5, legend_sw; end
    % axis square
end
fontsize(gcf, 14, "points")
saveas(gcf, "figs\figs_clean_all\rest_run_mean_roi_indmice.png")
saveas(gcf, "figs\figs_clean_all\rest_run_mean_roi_indmice.fig")