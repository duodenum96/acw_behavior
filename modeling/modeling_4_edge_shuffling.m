%% Shuffling of edges
clear, clc, close all
cd C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\modeling_3
addpath(genpath('C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\modeling\external\'))
addpath C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\

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
W(logical(eye(n))) = 0; % Important for shuffling

surrogates = 0:0.2:1; % Degree of shuffling
nsurro = length(surrogates);

nsims = 30;
ntp = 20;
acwdrs_rest = zeros(n, ntp, nsurro, nsims);
acwdrs_task = zeros(n, ntp, nsurro, nsims);
swp = zeros(nsurro, nsims);

parfor i = 1:nsims
    for j = 1:nsurro
        tic
        %% Set up W
        j_surr = surrogates(j);
        W_surro = shuffwei_partial(W, j_surr);
        % Average upper and lower triangles
        % Simulations
        W_i = W_surro;
        W_i(logical(eye(n))) = 1; % add recurrent connections
        %% run simulations

        x = wilsoncowan_RK2(tau, b, W_i, k, s, C, tspan);
        acwdrs_rest(:, :, j, i) = acw_windowed_loop_dr(x, 1/dt, 10);
    
        x = wilsoncowan_RK2(tau, b, W_i, k, s_task, C, tspan);
        acwdrs_task(:, :, j, i) = acw_windowed_loop_dr(x, 1/dt, 10);
        toc
        fprintf("\ni = %d, j = %d\n", i, j)
    end
end

save("cellreports_review\dtinetwork_edgeshuffling.mat", "acwdrs_rest", "acwdrs_task")
