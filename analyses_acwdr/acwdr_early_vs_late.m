clear, clc, close all
cd C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw
% cd E:\mice_regularacw
load("acwdr_ca.mat")
sessions = string(sessions);
mousenames = string(mousenames);

labels = ["onset", "offset", "initial rest", "sustained rest", "locomotion"];
labels_realdeal = strings(size(labels_onehot, 1), 1);
for i = 1:length(labels)
    indx = labels_onehot == i;
    labels_realdeal(indx) = labels(i);
end


susrest_acws = acwmat(labels_realdeal == "sustained rest");
loco_acws = acwmat(labels_realdeal == "locomotion");
susrest_sessions = sessions(labels_realdeal == "sustained rest");
loco_sessions = sessions(labels_realdeal == "locomotion");
susrest_mousenames = mousenames(labels_realdeal == "sustained rest");
loco_mousenames = mousenames(labels_realdeal == "locomotion");

mice = unique(mousenames);
nmice = length(mice);
susrest_initial = cell(1, nmice);
susrest_later = cell(1, nmice);
loco_initial = cell(1, nmice);
loco_later = cell(1, nmice);
for i = 1:nmice
    % Handle sustained rest
    susrest_initial{i} = [];
    susrest_later{i} = [];

    i_mice = mice(i);
    i_acw = susrest_acws(susrest_mousenames == i_mice);
    
    i_sessions = susrest_sessions(susrest_mousenames == i_mice);
    uniqses = unique(i_sessions);
    for j = 1:length(uniqses)
        j_sessions = i_sessions(i_sessions == uniqses(j));
        j_acw = i_acw(i_sessions == uniqses(j));
        if length(j_sessions) > 1
            susrest_initial{i} = [susrest_initial{i}; j_acw(1)];
            susrest_later{i} = [susrest_later{i}; j_acw(end)];
        end
    end
    
    % Now do locomotion
    loco_initial{i} = [];
    loco_later{i} = [];
    i_acw = susrest_acws(loco_mousenames == i_mice);
    
    i_sessions = loco_sessions(loco_mousenames == i_mice);
    uniqses = unique(i_sessions);
    for j = 1:length(uniqses)
        j_sessions = i_sessions(i_sessions == uniqses(j));
        j_acw = i_acw(i_sessions == uniqses(j));
        if length(j_sessions) > 1
            loco_initial{i} = [loco_initial{i}; j_acw(1)];
            loco_later{i} = [loco_later{i}; j_acw(end)];
        end
    end
end

% Do mice vector
susrest_init_micelengths = zeros(1, nmice);
susrest_late_micelengths = zeros(1, nmice);
loco_init_micelengths = zeros(1, nmice);
loco_late_micelengths = zeros(1, nmice);
for i = 1:nmice
    susrest_init_micelengths(i) = length(susrest_initial{i});
    susrest_late_micelengths(i) = length(susrest_later{i});
    loco_init_micelengths(i) = length(loco_initial{i});
    loco_late_micelengths(i) = length(loco_later{i});
end

susrest_micevector_init = micelengths2micevec(susrest_init_micelengths);
susrest_micevector_late = micelengths2micevec(susrest_late_micelengths);
loco_micevector_init = micelengths2micevec(loco_init_micelengths);
loco_micevector_late = micelengths2micevec(loco_late_micelengths);

susrest_initial_vec = cat(1, susrest_initial{:});
susrest_late_vec = cat(1, susrest_later{:});
loco_initial_vec = cat(1, loco_initial{:});
loco_late_vec = cat(1, loco_later{:});

%% Compare susrest
y = [susrest_initial_vec; susrest_late_vec];
x = [ones(length(susrest_initial_vec), 1); ones(length(susrest_late_vec), 1)*2];
mice_code = [susrest_micevector_init; susrest_micevector_late];

colors = slanCM('Set1', 5);
colmat = zeros(length(x), 3);
for i = 1:length(x)
    i_mice = mice_code(i);
    colmat(i, :) = colors(i_mice, :);
end

[~, ~, stats] = ranksum(susrest_initial_vec, susrest_late_vec);
[pval,r,U,rsum,h] = ranksum_effect_size(susrest_initial_vec, susrest_initial_vec);
stattext = get_stattext_wilcox(stats.zval, pval, r, 'horizontal');
swarmchart(x, y, [], colmat, 'filled')
xticks([1 2])
xticklabels(["Initial", "Late"])
sigstar_duodenal({[1 2]}, pval,stattext, 0)
ylabel("timescale (\tau)")
hold on
line([0.75 1.25], [median(y(x==1)) median(y(x==1))], 'color', 'k', 'LineWidth', 1.5)
line([1.75 2.25], [median(y(x==2)) median(y(x==2))], 'color', 'k', 'LineWidth', 1.5)
set(gca, 'YGrid', 'on')
legend_sw(colors)
title("Sustained Rest")
fontsize(gca, 14, "points")
% saveas(gcf, 'figs_acwdr\additional_analyses\early_vs_late_susrest.jpg')
exportgraphics(gcf, 'figs_acwdr\additional_analyses\early_vs_late_susrest.jpg', ...
    'Resolution', 600)
%% Compare loco
y = [loco_initial_vec; loco_late_vec];
x = [ones(length(loco_initial_vec), 1); ones(length(loco_late_vec), 1)*2];
mice_code = [loco_micevector_init; loco_micevector_late];

colors = slanCM('Set1', 5);
colmat = zeros(length(x), 3);
for i = 1:length(x)
    i_mice = mice_code(i);
    colmat(i, :) = colors(i_mice, :);
end

[~, ~, stats] = ranksum(loco_initial_vec, loco_late_vec);
[pval,r,U,rsum,h] = ranksum_effect_size(loco_initial_vec, loco_initial_vec);
stattext = get_stattext_wilcox(stats.zval, pval, r, 'horizontal');
swarmchart(x, y, [], colmat, 'filled')
xticks([1 2])
xticklabels(["Initial", "Late"])
sigstar_duodenal({[1 2]}, pval,stattext, 0)
ylabel("timescale (\tau)")
hold on
line([0.75 1.25], [median(y(x==1)) median(y(x==1))], 'color', 'k', 'LineWidth', 1.5)
line([1.75 2.25], [median(y(x==2)) median(y(x==2))], 'color', 'k', 'LineWidth', 1.5)
set(gca, 'YGrid', 'on')
legend_sw(colors)
title("Locomotion")
fontsize(gca, 14, "points")
exportgraphics(gcf, 'figs_acwdr\additional_analyses\early_vs_late_locomotion.jpg', ...
    'Resolution', 600)
