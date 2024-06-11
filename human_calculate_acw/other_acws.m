%% Organize other data
clear, clc, close all
addpath('C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\daviddata\task')
addpath('C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\')
cd C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\daviddata\task
files = dir('*other*_reordered.mat');
filenames = cat(1, files.name);
subjnames = filenames(:, 1:3);

filenames = string(filenames);
subjnames = string(subjnames);

n_subjs = length(subjnames);
n_chans = 64;
%% 10 second windows
windowsize = 10;

acw0s_other = cell(1, n_subjs);
acw50s_other = cell(1, n_subjs);

for i = 1:n_subjs
    tic
    eegdata = load(filenames(i));
    data = double(eegdata.EEG.data);
    fs = eegdata.EEG.srate;
    [acw0s_other{i}, acw50s_other{i}] = acw_windowed_loop(data, fs, windowsize);
    disp(i + " / " + n_subjs)
    toc
end

save("other_acws.mat", "acw50s_other", "acw0s_other", "subjnames")