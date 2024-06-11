%% Visualize human ML
clear, clc, close all
cd C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw
addpath(genpath('C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\ML_matlab_acw0'))
addpath(genpath('C:\Users\duodenum\Desktop\brain_stuff\misc'))

load('C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\ML_matlab_acw0\humans\human_mlresults.mat')

colors = [	"#0072BD" 	"#D95319" 	"#7E2F8E"];
labels = ["rest", "self", "other"];
%% Interpolate x axis values, match to each other
% roc_x: 10x5 cell (10 cvs, 5 labels)

xaxis = 0:0.01:1; % common x axis

% all_xs = zeros(length(xaxis), 10, 5);
all_ys = zeros(length(xaxis), 10, 5);
for i = 1:10
    for j = 1:3
        % roc_x = stufftoplot.roc_x{i,j};
        roc_x = roc_x_1vall{i,j};
        rocx_smooth = rocxsmooth(roc_x);
        roc_y = interp1(rocx_smooth, roc_y_1vall{i,j}, xaxis);
        all_ys(:, i, j) = roc_y;
    end
end
%% Plot average ROC curve +- shading as error bar
mean_auc = mean(cell2mat(stufftoplot.auc));
sd_auc = std(cell2mat(stufftoplot.auc));
% names = labels + " (M: " + mean_auc + ", STD: " + sd_auc + ")";
names = strings(3, 1);
for i = 1:3
    names(i) = sprintf("%s AUC: M = %.3f, STD =  %.3f", ...
        labels(i), mean_auc(i), sd_auc(i));
end

for i = 1:3
    shadedErrorBar(xaxis, all_ys(:, :, i)', {@mean,@std},...
        'lineprops', {'-', 'Color', colors(i), 'LineWidth', 3, ...
        'DisplayName', names(i)}, 'patchSaturation',0.15)
    hold on
end
xlim([-0.02 1.02])
ylim([0 1.05])
grid on
legend('Box', 'off')
xlabel('False Positive Rate')
ylabel('True Positive Rate')
fontsize(gcf, 14, "points")
% saveas(gcf, 'ML_matlab_acw0\mean_human_roc.jpg')
% saveas(gcf, 'ML_matlab_acw0\mean_human_roc.fig')
exportgraphics(gcf, 'ML_matlab_acw0\mean_human_roc.jpg', ...
    'Resolution', 600)
%% Plot individual ROC curves
close all
figure(Position=[4.2000 338 1.5296e+03 420]);
% tiledlayout(2, 5)
for i = 1:10
    subplot(2,5,i)
    for j = 1:3
        roc_x = stufftoplot.roc_x{i,j};
        roc_y = stufftoplot.roc_y{i,j};
        i_auc = stufftoplot.auc{i, j};
        name = sprintf("%s AUC: %.3f", labels(j), i_auc);
        plot(roc_x, roc_y, 'Color', ...
            colors(j), 'DisplayName', name)
        hold on
    end
    xlim([-0.02 1.02])
    ylim([0 1.05])
    xticks(0:0.2:1)
    grid on
    legend('Box', 'off', 'Location','southeast')
    axis square
    if i > 5, xlabel('False Positive Rate'), end
    if any(i == [1 6]), ylabel('True Positive Rate'), end
end
sgtitle("Folds")
% saveas(gcf, 'ML_matlab_acw0\mean_human_roc_allfolds.jpg')
% saveas(gcf, 'ML_matlab_acw0\mean_human_roc_allfolds.fig')
exportgraphics(gcf, 'ML_matlab_acw0\mean_human_roc_allfolds.jpg', ...
    'Resolution', 600)
%% Plot confusion matrices
% All data
close all
all_confusion = sum(cat(3, stufftoplot.confusion{:}), 3);
% yaxis: true class (rows); xaxis: predicted class (columns)
confusionchart(all_confusion, labels, ...
    'RowSummary','row-normalized')
fontsize(gcf, 14, "points")
% saveas(gcf, 'ML_matlab_acw0\confusionmatrix_all_human.fig')
% saveas(gcf, 'ML_matlab_acw0\confusionmatrix_all_human.jpg')
exportgraphics(gcf, 'ML_matlab_acw0\confusionmatrix_all_human.jpg', ...
    'Resolution', 600)

% Averaged
% close all
% mean_confusion = mean(cat(3, stufftoplot.confusion{:}), 3);
% sd_confusion = std(cat(3, stufftoplot.confusion{:}), [], 3);
% stringmat = cell(3);
% for i = 1:3
%     for j = 1:3
%         stringmat{i, j} = sprintf('M: %.3f\nSTD: %.3f', ...
%             mean_confusion(i, j), sd_confusion(i, j));
%     end
% end
% heatmap(mean_confusion, labels, labels, stringmat, ...
%     'textcolor', 'r')
% xlabel('Predicted Class')
% ylabel('True Class')
% saveas(gcf, 'ML_matlab\confusionmatrix_mean_human.fig')
% saveas(gcf, 'ML_matlab\confusionmatrix_mean_human.jpg')
%% Plot all confusion matrices
% tiledlayout(2,5)
close all
for i = 1:10
    subplot(2, 5, i)
    confusionchart(stufftoplot.confusion{i}, labels)
end
sgtitle("Folds")
% saveas(gcf, 'ML_matlab_acw0\confusionmatrix_allfolds_human.fig')
% saveas(gcf, 'ML_matlab_acw0\confusionmatrix_allfolds_human.jpg')
exportgraphics(gcf, 'ML_matlab_acw0\confusionmatrix_allfolds_human.jpg', ...
    'Resolution', 600)

%% Helper function
function rocx_mod = rocxsmooth(roc_x)
    rocx_mod = roc_x;
    for i = 2:length(rocx_mod)
        if rocx_mod(i) <= rocx_mod(i-1)
            rocx_mod(i) = rocx_mod(i) + abs(rocx_mod(i) - rocx_mod(i-1)) + eps;
        end
    end
end
