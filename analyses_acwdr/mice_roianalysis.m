clear, clc, close all
addpath(genpath('C:\Users\duodenum\Desktop\brain_stuff\misc'))
cd C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw
load("acwdr_ca.mat")
rois = load("revision\mice_rois.mat");

labels = ["onset", "offset", "initial rest", "sustained rest", "locomotion"];
labels_realdeal = strings(size(labels_onehot, 1), 1);
for i = 1:length(labels)
    indx = labels_onehot == i;
    labels_realdeal(indx) = labels(i);
end

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

js = [1 5 2 3 4]; % onset, run, offset, initrest, susrest

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
% Just copy - pasted a bunch of lines from the pipeline
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
acw_var_x_trials_running = zeros(nroi, nmice);
acw_mean_x_trials_running = zeros(nroi, nmice);
for i = 1:nmice
    mouse = string(mice{i});
    mouseindx = restmice == mouse;
    acws = restacws(mouseindx, :);
    acw_var_x_trials(:, i) = std(acws);
    acw_mean_x_trials(:, i) = mean(acws);

    mouseindx = runmice == mouse;
    acws = runacws(mouseindx, :);
    acw_var_x_trials_running(:, i) = std(acws);
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
    ntrials_running(i) = length(acw_mean_x_rois_running{i});
end
% old name: mice_onehot
mice_onehot_x_roi = [ones(ntrials(1), 1); 2*ones(ntrials(2), 1); ...
    3*ones(ntrials(3), 1); 4*ones(ntrials(4), 1); 5*ones(ntrials(5), 1)];

mice_onehot_x_roi_run = [ones(ntrials_running(1), 1); 2*ones(ntrials_running(2), 1); ...
    3*ones(ntrials_running(3), 1); 4*ones(ntrials_running(4), 1); 5*ones(ntrials_running(5), 1)];

 % In the end, we have acw_var_x_rois, acw_mean_x_rois, acw_var_x_trials,
 % acw_mean_x_trial and the *_running variants of these

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

% save('restvariability_acwdr.mat', "calculations")

calculations = [];
calculations.acw_mean_x_rois_running = acw_mean_x_rois_running;
calculations.acw_mean_x_trials_running = acw_mean_x_trials_running;
calculations.acw_mean_x_rois_running_col = cat(1, acw_mean_x_rois_running{:});
calculations.mice_onehot_x_roi_run = mice_onehot_x_roi_run;

% save('runvariability_acwdr.mat', "calculations")

%% Calculate prc change of mean and var
acw_mean_prc_x_trl = ((acw_mean_x_trials_running - acw_mean_x_trials) ./ ...
    acw_mean_x_trials) .* 100;

%% ACTUAL ANALYSIS
roigroups = fieldnames(rois);
n_roigroups = length(roigroups);
name_roigroups = ["frontal", "motor", ...
    "somatosensory1", "somatosensory2", "somatosensory3", ...
    "visual"];
name_roigroups2 = [name_roigroups + "_L", name_roigroups + "_R"];
names_nice = ["Left Frontal", "Left Motor", ...
    "Left Somatosensory 1", "Left Somatosensory 2", "Left Somatosensory 3", ...
    "Left Visual", ...
    "Right Frontal", "Right Motor", ...
    "Right Somatosensory 1", "Right Somatosensory 2", "Right Somatosensory 3", ...
    "Right Visual"];


%% Rest - Task Difference
%% Get MORE colors (because why not)
colors_sd = slanCM('Set3', 5);

ncol = size(mice_onehot_x_trial, 1);
colmat_sd = zeros(ncol, 3);
for i = 1:nmice
    colmat_sd(mice_onehot_x_trial == i, :) = repmat(colors_sd(i, :), nroi, 1);
end
load('acwdr_ca.mat', 'sessions')
runtrl = labels_names == "running"; 
resttrl = labels_names == "sustained rest";

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
rest_stds = zeros(n_trials, nroi);
mousenames_averages_rest =  zeros(n_trials, 1);

for i =1:n_trials
    indx = restsessions == uniqtrials(i);
    tmp = mousenames_rest_code(indx, 1);
    mousenames_averages_rest(i) = tmp(1);
    rest_averages(i,:) = mean(acw_rest(indx, :), 1); % For each ROI, average across trials
    % rest_stds(i,:) = std(acw_rest(indx, :), [], 1); % For each ROI, average across trials
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
    % rest_stds(i,:) = std(acw_run(indx, :), [], 1); % For each ROI, average across trials
end
nrois = 92;
%% Multiple Comparisons Correction
pvals = zeros(n_roigroups, 1);
for j = 1:n_roigroups
    i_roigroup = roigroups(j);
    i_rois = rois.(i_roigroup{1});
    i_acw_var_x_trials = acw_var_x_trials(i_rois, :);
    i_acw_mean_prc_x_trl = acw_mean_prc_x_trl(i_rois, :);
    i_mice_onehot_x_trial = repmat(1:5, length(i_rois), 1);
    
    i_acw_rest = acw_rest(:, i_rois);
    i_acw_run = acw_run(:, i_rois);

    restacw_ave = zeros(length(i_rois), nmice);
    runacw_ave = zeros(length(i_rois), nmice);
    for i = 1:nmice
        i_restmice = mousenames_rest_code(:, 1) == i;
        i_runmice = mousenames_run_code(:, 1) == i;
        restacw_ave(:, i) = mean(i_acw_rest(i_restmice, :), 1);
        runacw_ave(:, i) = mean(i_acw_run(i_runmice, :), 1);
    end

    [~, ~, stats] = ranksum(restacw_ave(:), runacw_ave(:));
    [pvals(j),r,U,rsum,h] = ranksum_effect_size(restacw_ave(:), runacw_ave(:));
end
pvals_correct = bonf_holm(pvals);
%% PLOT
for j = 1:n_roigroups
    close all
    pval = pvals_correct(j);

    i_roigroup = roigroups(j);
    i_rois = rois.(i_roigroup{1});
    i_acw_var_x_trials = acw_var_x_trials(i_rois, :);
    i_acw_mean_prc_x_trl = acw_mean_prc_x_trl(i_rois, :);
    i_mice_onehot_x_trial = repmat(1:5, length(i_rois), 1);
    
    i_acw_rest = acw_rest(:, i_rois);
    i_acw_run = acw_run(:, i_rois);

    restacw_ave = zeros(length(i_rois), nmice);
    runacw_ave = zeros(length(i_rois), nmice);
    for i = 1:nmice
        i_restmice = mousenames_rest_code(:, 1) == i;
        i_runmice = mousenames_run_code(:, 1) == i;
        restacw_ave(:, i) = mean(i_acw_rest(i_restmice, :), 1);
        runacw_ave(:, i) = mean(i_acw_run(i_runmice, :), 1);
    end
    mice_code = repmat(1:5, length(i_rois), 1);
    
    y = [restacw_ave(:); runacw_ave(:)];
    x = [ones(length(restacw_ave(:)), 1); 2 * ones(length(restacw_ave(:)), 1)];
    mice_code2 = [mice_code(:); mice_code(:)];
    
    [~, ~, stats] = ranksum(restacw_ave(:), runacw_ave(:));
    [~,r,U,rsum,h] = ranksum_effect_size(restacw_ave(:), runacw_ave(:));
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
    sigstar_duodenal({[1 2]}, pval,stattext, 0, 0.05)
    ylabel("Mean of timescale (\tau) across windows")
    hold on
    line([0.75 1.25], [median(y(x==1)) median(y(x==1))], 'color', 'k', 'LineWidth', 1.5)
    line([1.75 2.25], [median(y(x==2)) median(y(x==2))], 'color', 'k', 'LineWidth', 1.5)
    set(gca, 'YGrid', 'on')
    legend_sw(colors, 'NorthWest')
    % title(mice(i))
    % xtickangle(15)
    
    fontsize(gcf, 16, "points")
    % saveas(gcf, "figs_acw50\resttaskdiff\mean_x_trl_allmice.jpg")
    % saveas(gcf, "figs_acw50\resttaskdiff\mean_x_trl_allmice.fig")
    title(names_nice(j))
    pause
    exportgraphics(gcf, ...
        sprintf("figs_acwdr\\resttaskdiff\\rois\\mean_x_trl_allmice_%s.jpg", name_roigroups2(j)), ...
        'Resolution',600)
    close all
end

