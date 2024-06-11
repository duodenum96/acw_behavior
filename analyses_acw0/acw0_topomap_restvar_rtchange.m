%% Correlation between ROI ordering and prc change
% c) Correlation between ROI ordering and ACW-0
% percent change from rest to running
%% Load stuff
clear, clc, close all
addpath('C:\Users\duodenum\Desktop\brain_stuff\misc\hex_and_rgb')
addpath C:\Users\duodenum\Desktop\brain_stuff\misc\slanCM

cd C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw
load topo_acw0.mat
load('yaxis_info.mat', 'order')
calculations = load('restvariability_acwdr.mat', "calculations");
acw_var_x_trials = calculations.calculations.acw_var_x_trials;

load('acw0_ca.mat', 'sessions');
sessions = string(sessions);

labels = ["onset", "offset", "initial rest", "sustained rest", "running"];

labels_names = strings(size(labels_onehot, 1), 1);
for i = 1:length(labels)
    indx = labels_onehot == i;
    labels_names(indx) = labels(i);
end

dict = table(labels_names, acwmat, string(mousenames), sessions, ...
    'VariableNames', {'labels', 'acw', 'mice', 'sessions'});

susrest = dict.labels == "sustained rest";
run = dict.labels == "running";
sesh = dict.sessions;
%% Dependent variables
rest = dict(susrest, :);
restacws = rest.acw;
restmice = rest.mice;
restses = sesh(susrest);

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

    acw_prcchange_x_trials(:, i) = (i_runacws_ave - i_restacws_ave) ./ i_restacws_ave;
end
%% Get colors
colors = [102,194,165
    252,141,98
    141,160,203
    231,138,195
    166,216,84] ./ 255;

%% Plot topographies of rest task change and rest variability
%% Import one topo for each mouse

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

%% Assign topographies for each state

% nlabels = length(labels);
acwtopos = cell(1, nmice);
n_thingstoplot = 2; % 1: variability, 2: prc change
for i = 1:nmice
    mouse = mice{i};
    % mouse_indx = string(mousenames) == string(mouse);
    acwtopos{i} = cell(1, n_thingstoplot);
    for j = 1:n_thingstoplot
        acwtopos{i}{j} = zeros(roidims);
        % label_indx = labels_onehot == j;
        % acwtopo = mean(acwmat(mouse_indx & label_indx, :));
        for k = 1:nroi
            roi_indx = topo{i} == k;
            if j == 1
                acwtopos{i}{j}(roi_indx) = acw_var_x_trials(k, i);
            elseif j == 2
                acwtopos{i}{j}(roi_indx) = acw_prcchange_x_trials(k, i);
            end
        end
    end
end
% Set nonrois to nan
for i = 1:nmice
    for j = 1:n_thingstoplot
        indx = acwtopos{i}{j} == 0;
        acwtopos{i}{j}(indx) = nan;
    end
end
%% Nice topoplots
labels = {"Rest variability of ACW-0 (s)", ["% ACW-0 change from", "rest to locomotion"]};
close all
figure(Position= [488 100 560 658]);
tiledlayout(nmice, n_thingstoplot, "TileSpacing","tight", "Padding", "tight")
for i = 1:nmice
    for j  = 1:n_thingstoplot
        ax = nexttile;
        if j == 1
            h = imagesc(acwtopos{i}{j});
            colormap(ax, 'summer')
            colorbar
        elseif j == 2
            h = imagesc(acwtopos{i}{j} .* 100);
            clims = [min([acwtopos{i}{j} .* 100], [],  "all"), max([acwtopos{i}{j} .* 100], [],  "all")];
            climmax = max(abs(clims));
            clims_symmetric = [-climmax, climmax];
            clim(clims_symmetric)
    
            cb = colorbar;
            % ylabel(cb,'% Change of STD (Loco - Rest)')
            set(cb, 'Limits', clims)
            colormap(ax, slanCM('vik'))
        end
        set(h, 'AlphaData', 1 - isnan(acwtopos{i}{j})) % get rid of background
        % clim([min([acwtopos{i}{:}], [],  "all"), max([acwtopos{i}{:}], [],  "all")])
        axis square
        axis off
        if i == 1, title(labels{j}), end
        if j == 1, ylabel(mice(i)), end
        % to check
        % title(mice(i) + labels(j))
    end
end
fontsize(gcf, 13, "points")
% saveas(gcf, 'figs_acw50\topomaps_variability_prchange.jpg')
% saveas(gcf, 'figs_acw50\topomaps_variability_prchange.fig')
exportgraphics(gcf, 'figs_acw0\topomaps_variability_prchange.jpg', ...
    'Resolution', 600)
close
save('restvar_varchange_data.mat')