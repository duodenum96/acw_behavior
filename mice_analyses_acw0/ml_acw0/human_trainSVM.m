clear, clc, close all
cd C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\daviddata\
addpath C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\ML_matlab_acw0\
load("acwdata_acw0.mat")

[trainedClassifier, validationAccuracy, stufftoplot, roc_x_1vall, roc_y_1vall] = ...
    trainSVM_human(all_acws, labels);

filename = ...
    'C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\ML_matlab_acwdr\humans\human_mlresults.mat';

save(filename, ...
    'trainedClassifier', 'validationAccuracy', "stufftoplot", ...
    "roc_x_1vall", "roc_y_1vall")
disp('Done')

