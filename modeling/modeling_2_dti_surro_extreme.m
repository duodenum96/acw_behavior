%% Modeling: Clustered network
clear, clc, close all
cd C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\modeling_3\
addpath C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\modeling\external\
addpath C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\
%% Parameters
b = -3; 
tau = 0.1;
C = 1;
s = 0; % rest
s_task = 1;
k = 1;
tspan = 300;
dt = 10 / 1000;

W = load('C:\Users\duodenum\Desktop\brain_stuff\modelingresources\averageConnectivity_Fpt.mat');
W = 10.^(W.Fpt);
n = length(W);
W(~logical(eye(n))) = (2*360)*(W(~logical(eye(n))) ./ sum(W(~logical(eye(n))), 'all'));

surrogates = 0:0.5:5;
nsurro = length(surrogates);

nsims = 30;
ntp = 20;
acwdrs_rest = zeros(n, ntp, nsurro, nsims);
acwdrs_task = zeros(n, ntp, nsurro, nsims);

mean_rest = zeros(n, nsurro, nsims);
mean_task = zeros(n, nsurro, nsims);
var_rest = zeros(n, nsurro, nsims);
var_task = zeros(n, nsurro, nsims);

parfor i = 1:nsims
    for j = 1:nsurro
        tic
        %% Set up W
        W_i = W;
        W_i(logical(eye(n))) = surrogates(j); % this is important before shuffling
        %% run simulations

        x = wilsoncowan_RK2(tau, b, W_i, k, s, C, tspan);
        acwdrs_rest(:, :, j, i) = ...
            acw_windowed_loop_dr(x, 1/dt, 10);
        mean_rest(:, j, i) = mean(x, 2);
        var_rest(:, j, i) = std(x, [], 2);
    
        x = wilsoncowan_RK2(tau, b, W_i, k, s_task, C, tspan);
        acwdrs_task(:, :, j, i) = ...
            acw_windowed_loop_dr(x, 1/dt, 10);
        mean_task(:, j, i) = mean(x, 2);
        var_task(:, j, i) = std(x, [], 2);
        toc
        fprintf("\ni = %d, j = %d\n", i, j)
    end
end

save("cellreports_review\dtinetwork_recurrent_extreme.mat", ...
    "acwdrs_rest", "acwdrs_task", "var_task", ...
    "mean_task", "var_rest", "mean_rest")


