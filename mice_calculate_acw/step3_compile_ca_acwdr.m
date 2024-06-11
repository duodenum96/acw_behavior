%% Compile all ACWs and put in matrices
clear
clc
close all

addpath('/BICNAS2/ycatal/mice_dynacw/')
addpath('/BICNAS2/ycatal/mice_dynacw/scripts/acws_per_state/auxiliary')
rootDIR = '/BICNAS2/ycatal/mice_dynacw/';
dataDIR = '/BICNAS2/ycatal/mice_data/WFOM_extracted_signals/'; % ROI data
boutsDIR = '/BICNAS2/ycatal/mice_data/beh_ext/'; 
acwDIR = '/BICNAS2/ycatal/mice_data/acw_states/'; % load all blocks
stateDIR =  '/BICNAS2/ycatal/mice_data/acw_states/states'; 

mice = {'cm124','cm125','cm126','cm127','cm128'};
states = {'offset', 'onset', 'post_offset', 'running'};

all_acw = [];
all_acw.onset = [];
all_acw.offset = [];
all_acw.initrest = [];
all_acw.susrest = [];
all_acw.running = [];

all_acw.mousenames_onset = [];
all_acw.mousenames_offset = [];
all_acw.mousenames_initrest = [];
all_acw.mousenames_susrest = [];
all_acw.mousenames_running = [];

all_acw.session_onset = [];
all_acw.session_offset = [];
all_acw.session_initrest = [];
all_acw.session_susrest = [];
all_acw.session_running = [];

for m = 1:length(mice)
    mouse = mice{m};
    for i = 1:length(states)
        s = states{i};
        data = load(fullfile(acwDIR, s, mouse));
        data = data.runningBlocks;
        nruns = length(data);

        if string(s) == string('post_offset')
            for j = 1:nruns
                all_acw.initrest = [all_acw.initrest, data(j).initrest_acw.ca_acwdr];
                all_acw.susrest = [all_acw.susrest, data(j).susrest_acw.ca_acwdr];

                all_acw.mousenames_initrest = [all_acw.mousenames_initrest; data(j).mouse_name];
                all_acw.mousenames_susrest = [all_acw.mousenames_susrest; data(j).mouse_name];

                all_acw.session_initrest = [all_acw.session_initrest; data(j).session];
                all_acw.session_susrest = [all_acw.session_susrest; data(j).session];
            end
        else
            for j = 1:nruns
                all_acw.(s) = [all_acw.(s), data(j).([s, '_acw']).ca_acwdr];
                all_acw.(['mousenames_', s]) = [all_acw.(['mousenames_', s]); data(j).mouse_name];
                all_acw.(['session_', s]) = [all_acw.(['session_', s]); data(j).session];
            end
        end
    end
end

labels = {'onset', 'offset', 'initrest', 'susrest', 'running'};
labels_onehot = [ones(size(all_acw.onset, 2), 1); ...
                 2 * ones(size(all_acw.offset, 2), 1); ...
                 3 * ones(size(all_acw.initrest, 2), 1); ...
                 4 * ones(size(all_acw.susrest, 2), 1); ...
                 5 * ones(size(all_acw.running, 2), 1)];

acwmat = [all_acw.onset'; all_acw.offset'; all_acw.initrest'; all_acw.susrest'; all_acw.running'];
mousenames = [all_acw.mousenames_onset; all_acw.mousenames_offset; all_acw.mousenames_initrest; ...
              all_acw.mousenames_susrest; all_acw.mousenames_running];
sessions = [all_acw.session_onset; all_acw.session_offset; all_acw.session_initrest; ...
    all_acw.session_susrest; all_acw.session_running];

save(fullfile(acwDIR, 'states', 'acwdr_ca'), 'acwmat', 'labels', ...
    'labels_onehot', 'all_acw', 'mousenames', 'sessions')

    