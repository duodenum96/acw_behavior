%% Load stuff
clear, clc, close all

cd C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw
load topo_acwdr.mat
load('yaxis_info.mat', 'order')

load('acwdr_ca.mat', 'sessions');
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
% Independent variables:
%   Resting state ACW variability across trials
%
% Dependent variables:
%   Rest - Running ACW percent change 


susrest = dict.labels == "sustained rest";
run = dict.labels == "running";
sesh = dict.sessions;

%% independet: Rest ACW variability across trials
rest = dict(susrest, :);
restacws = rest.acw;
restmice = rest.mice;
restses = sesh(susrest);

% For every mouse, calculate variability of every ROI across trials
acw_var_x_trials = zeros(nroi, nmice);
acw_mean_x_trials = zeros(nroi, nmice);
for i = 1:nmice
    mouse = string(mice{i});
    mouseindx = restmice == mouse;
    acws = restacws(mouseindx, :);
    acw_var_x_trials(:, i) = std(acws);
    acw_mean_x_trials(:, i) = mean(acws);
end
mice_onehot_rois = [ones(nroi, 1); 2*ones(nroi, 1); 3*ones(nroi, 1); 4*ones(nroi, 1); 5*ones(nroi, 1)];

% Same for across ROIs
acw_mean_x_rois = cell(nmice, 1);
ntrials = zeros(nmice, 1);
for i = 1:nmice
    mouse = string(mice{i});
    mouseindx = restmice == mouse;
    acws = restacws(mouseindx, :);
    acw_mean_x_rois{i} = mean(acws, 2);
    ntrials(i) = length(acw_mean_x_rois{i});
end
mice_onehot = [ones(ntrials(1), 1); 2*ones(ntrials(2), 1); ...
    3*ones(ntrials(3), 1); 4*ones(ntrials(4), 1); 5*ones(ntrials(5), 1)];

%% Get colors
colors = [102,194,165
          252,141,98
          141,160,203
          231,138,195
          166,216,84] ./ 255;
ncol = size(mice_onehot_rois, 1);
colmat = zeros(ncol, 3);
for i = 1:nmice
    colmat(mice_onehot_rois == i, :) = repmat(colors(i, :), nroi, 1);
end
%% save

calculations = [];
calculations.acw_var_x_trials = acw_var_x_trials;
calculations.acw_var_x_trials_mice = mice_onehot_rois;
calculations.acw_var_x_trials_col = acw_var_x_trials(:);
calculations.acw_mean_x_rois = acw_mean_x_rois;

save('restvariability_acwdr.mat', "calculations")
%% Dependent variables
rundata = dict(run, :);
runacws = rundata.acw;
runmice = rundata.mice;
runses = sesh(run);

% For every ROI, calculate average rest acw, average run acw and 
% calculate prc change
acw_prcchange_x_trials = zeros(nroi, nmice);
for i = 1:nmice
    mouse = string(mice{i});
    mouseindx = runmice == mouse;
    i_runacws_ave = mean(runacws(mouseindx, :));

    mouseindx = restmice == mouse;
    i_restacws_ave = mean(restacws(mouseindx, :));

    acw_prcchange_x_trials(:, i) = ((i_runacws_ave - i_restacws_ave) ./ i_restacws_ave) .* 100;
end
%% Save percent change
calculations = [];
calculations.acw_prcchange_x_trials = acw_prcchange_x_trials;
calculations.mice_onehot_rois = mice_onehot_rois;

save('rest_run_prcchange_acwdr.mat', "calculations")
%% Calculate running ACW variability

% For every mouse, calculate variability of every ROI across trials during
% running

% Same for across ROIs
acw_mean_x_rois_running = cell(nmice, 1);
ntrials = zeros(nmice, 1);
for i = 1:nmice
    mouse = string(mice{i});
    mouseindx = runmice == mouse;
    acws = runacws(mouseindx, :);
    acw_mean_x_rois_running{i} = mean(acws, 2);
    ntrials(i) = length(acw_mean_x_rois_running{i});
end
mice_onehot = [ones(ntrials(1), 1); 2*ones(ntrials(2), 1); ...
    3*ones(ntrials(3), 1); 4*ones(ntrials(4), 1); 5*ones(ntrials(5), 1)];

%% save

calculations = [];
calculations.acw_var_x_rois_mice_running = mice_onehot;
calculations.acw_var_x_trials_mice_running = mice_onehot_rois;
calculations.acw_mean_x_rois_col_running = cat(1, acw_mean_x_rois_running{:});
calculations.acw_mean_x_rois_running = acw_mean_x_rois_running;

save('runvariability_acwdr.mat', "calculations")
