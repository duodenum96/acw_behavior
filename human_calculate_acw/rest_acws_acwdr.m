%% Organize Rest data
clear, clc, close all
addpath C:\Users\duodenum\Desktop\brain_stuff\misc
addpath('C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\daviddata\rest')
addpath('C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\')
addpath 'C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\daviddata'
cd C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\daviddata\rest
files = dir('*reordered.mat');
filenames = cat(1, files.name);
subjnames = filenames(:, 1:3);

filenames = string(filenames);
subjnames = string(subjnames);

n_subjs = length(subjnames);
n_chans = 64;
%% 10 second windows
windowsize = 10;

acwdrs_rest = cell(1, n_subjs);

for i = 1:n_subjs
    tic
    eegdata = load(filenames(i));
    data = double(eegdata.EEG.data);
    fs = eegdata.EEG.srate;
    [acwdrs_rest{i}] = acw_windowed_loop_dr(data, fs, windowsize);
    disp(i + " / " + n_subjs)
    toc
end

save("rest_acws_acwdr.mat", "acwdrs_rest", "subjnames")