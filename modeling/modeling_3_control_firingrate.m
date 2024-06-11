%% Modeling: Clustered network
clear, clc, close all
cd C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\modeling_3\
addpath C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\modeling\external\
addpath C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\
addpath 'C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\daviddata'
%% Parameters
b = -3; 
tau = 0.1;
C = 1;
k = 1;
tspan = 300;
dt = 10 / 1000;
chunk = 50;
burn = 20;

rest_target = 0.1; % Targeted firing rate
task_target = 0.6; % Targeted firing rate
s_rest_init = 0;   % Initial guess
s_task_init = 1;

W = load('C:\Users\duodenum\Desktop\brain_stuff\modelingresources\averageConnectivity_Fpt.mat');
W = 10.^(W.Fpt);
n = length(W);
W(~logical(eye(n))) = (2*360)*(W(~logical(eye(n))) ./ sum(W(~logical(eye(n))), 'all'));

surrogates = [0 0.5 1 1.5 2 2.5 3]; % TODO: CHANGE THIS
nsurro = length(surrogates);

nsims = 30; % TODO: CHANGE THIS
ntp = 20;
acwdrs_rest = zeros(n, ntp, nsurro, nsims);
acwdrs_task = zeros(n, ntp, nsurro, nsims);
mean_rest = zeros(n, nsurro, nsims);
mean_task = zeros(n, nsurro, nsims);
var_rest = zeros(n, nsurro, nsims);
var_task = zeros(n, nsurro, nsims);

optimal_s_rest = zeros(1, nsurro);
optimal_s_task = zeros(1, nsurro);
% Set the appropriate input values for different recurrent connections
for i = 1:nsurro
    W_i = W;
    W_i(logical(eye(n))) = surrogates(i); % this is important before shuffling
    optimal_s_rest(i) = wilsoncowan_RK2_graddescent(tau, b, W_i, k, ...
        s_rest_init, C, chunk, dt, burn, rest_target);
    optimal_s_task(i) = wilsoncowan_RK2_graddescent(tau, b, W_i, k, ...
        s_task_init, C, chunk, dt, burn, task_target);
end

parfor i = 1:nsims
    for j = 1:nsurro
        tic
        %% Set up W
        i_rec = surrogates(j);
        W_i = W;
        W_i(logical(eye(n))) = i_rec; % this is important before shuffling
        %% run simulations

        x = wilsoncowan_RK2(tau, b, W_i, k, optimal_s_rest(j), C, tspan);
        acwdrs_rest(:, :, j, i) = acw_windowed_loop_dr(x, 1/dt, 10);
        mean_rest(:, j, i) = mean(x, 2);
        var_rest(:, j, i) = std(x, [], 2);
 
    
        x = wilsoncowan_RK2(tau, b, W_i, k, optimal_s_task(j), C, tspan);
        acwdrs_task(:, :, j, i) = acw_windowed_loop_dr(x, 1/dt, 10);
        mean_task(:, j, i) = mean(x, 2);
        var_task(:, j, i) = std(x, [], 2);
        toc
        fprintf("\ni = %d, j = %d\n", i, j)
    end
end

save("cellreports_review\dtinetwork_control_fr.mat", "acwdrs_rest", "acwdrs_task", "var_task", ...
    "mean_task", "var_rest", "mean_rest")
save("cellreports_review\optimal_strengths.mat", "surrogates", "optimal_s_task", "optimal_s_rest")