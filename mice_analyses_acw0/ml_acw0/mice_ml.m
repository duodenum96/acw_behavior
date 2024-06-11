%% Train mice ML

clear, clc, close all
cd C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw
addpath 'C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\ML_matlab_acw0'
load("acw0_ca.mat")

labels = ["onset", "offset", "initial rest", "sustained rest", "locomotion"];
labels_realdeal = strings(size(labels_onehot, 1), 1);
for i = 1:length(labels)
    indx = labels_onehot == i;
    labels_realdeal(indx) = labels(i);
end

[trainedClassifier, validationAccuracy, stufftoplot, roc_x_1vall, roc_y_1vall] = ...
    trainSVM(acwmat, labels_realdeal);

save('ML_matlab_acw0\mice_mlresults.mat', ...
    'trainedClassifier', 'validationAccuracy', "stufftoplot", "roc_x_1vall", "roc_y_1vall")
disp('Done')