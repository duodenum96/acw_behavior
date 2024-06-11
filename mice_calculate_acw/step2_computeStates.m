% This script generates behavioral states for each mouse, similar to the 
% ones depicted in Figures 3B, 7A, S4, S6, and S7. The correlation maps are 
% calculated over a 10-second window, corresponding to 200 frames.
% 
% Author: Somayeh "Bahar" Shahsavarani
% Modified by Yasir
% email: bahar@huskers.unl.edu

%% initialize the directories
clear;clc;

%{
% rootDIR = '/home/ss6238/hillman_servers/';
% dataDIR = strcat(rootDIR,'enterprise/3/Bahar/organizing_data/data/');  % ROI data
% boutsDIR = strcat(rootDIR,'enterprise/3/Bahar/states_025/24-codes4share/Analysis/Sanity_Testing/runningBlocks/'); % running blocks created in step1
% acwDIR = strcat(rootDIR,'enterprise/3/Bahar/states_025/24-codes4share/Analysis/Sanity_Testing/states/runningBlocksWCorrMaps/');  % to save running blocks with correlation maps
% stateDIR = strcat(rootDIR,'enterprise/3/Bahar/states_025/24-codes4share/Analysis/Sanity_Testing/states/'); % to save behavioral states
%}

addpath('/BICNAS2/ycatal/mice_dynacw/')
addpath('/BICNAS2/ycatal/mice_dynacw/scripts/behavioral/auxiliary')
addpath('/BICNAS2/ycatal/mice_dynacw/scripts/acws_per_state/auxiliary')
rootDIR = '/BICNAS2/ycatal/mice_dynacw/';
dataDIR = '/BICNAS2/ycatal/mice_data/WFOM_extracted_signals/'; % ROI data
boutsDIR = '/BICNAS2/ycatal/mice_data/beh_ext/'; % running blocks created in step1
acwDIR = '/BICNAS2/ycatal/mice_data/acw_states/'; % to save running blocks with ACWs
stateDIR =  '/BICNAS2/ycatal/mice_data/acw_states/states'; % to save behavioral states


mice = {'cm124','cm125','cm126','cm127','cm128'};
fs = 20;

%% Compute correlation maps related to each state
%
% Correlation maps over each behavioral state are calculated for both neural
% and hemodynamic acitivty.
% The average of these maps is considered as the states.
% Five states are identified based on the mouse behavior: locomotion onset,
% locomotion (runnig), locomotion offset, resting states (post-offset) including
% initial rest (right after mouse stops running) and sustained state (when
% mouse has not been running for a while, i.e., 40 s), locomotion onset
%
% Within acwDIR, you need to have 4 folders called: onset, running,
% offset, and post_offset. The struct running blocks with new field related
% to the ACW maps for each state will be saved here.
%
% To be able to run this section, you need following auxiliary codes:
%       - running_acws
%       - offset_acws
%       - postoffset_acws
%       - onset_acws


delay = 30; % 1.5 s delay for hemodynamic data
for i = 1:length(mice)
    mousename = mice{i};
    disp('Processing: ')
    disp(mousename)
    %
    % load running bouts
    load(strcat(boutsDIR,mousename,'.mat'))
    
    %
    % compute running ACW maps
    running_acws(mousename,runningBlocks,acwDIR,dataDIR,delay,fs);
    
    % compute offset ACW maps
    offset_acws(mousename,runningBlocks,acwDIR,dataDIR,delay, fs);
    %
    % comput post-running ACW maps
    % 10-s window
    postoffset_acws(mousename,runningBlocks,acwDIR,dataDIR,delay, fs);
    %
    % compute onset ACW maps
    onset_acws(mousename,runningBlocks,acwDIR,dataDIR,delay, fs);
    %
    
    clear runningBlocks runningBlocks2 runningBlocks3 runningBlocks4 runningBlocks_final
end

disp('DONE')
