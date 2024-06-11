function runningBlocks = postoffset_acws(mousename,runningBlocks,acwDIR,dataDIR,delay, fs)
    % This function computes the correlation maps over the resting-state windows 
    % for each running block. The resting states are initial rest and sustained
    % rest. The initial rest is defined as 10 s rest immediately after
    % locomotion offset. The sustained rest is defined as 10 s rest starting
    % from 40 s after locomotion offset.
    % For each mouse, it returns a struct of running blocks with two new fields: 
    % correlations_neural and correlations_HbT, which contain correlation maps
    % over the resting-state winows for neural and hemo data.
    % 
    % Author: Somayeh "Bahar" Shahsavarani
    % email: bahar@huskers.unl.edu
    %
    % Modified by Yasir for ACW values
    
    % remove running bouts with less than 5 s duration
    duration = [runningBlocks.duration];
    durationIDX = duration < 100;
    runningBlocks(durationIDX) = [];
    
    % remove running bouts with post-running rest time less than 60 s
    time2nextrun = [runningBlocks.time2nextrun];
    time2nextrunIDX = time2nextrun < 1200;
    runningBlocks(time2nextrunIDX) = [];
    
    startFrame1 = 1; % the start frame affter offset for the initial rest (immediately after run ends)
    startFrame2 = 800; % the start frame affter offset for the sustained rest (40 s after run ends)
    
    % for FC maps of 10 seconds, set it to 100 for half a window
    nframes = 100;
    dt = 1 / fs;
    for run = 1:length(runningBlocks)
        sessionName = runningBlocks(run).session;
        
        fprintf('run %d ... \n', run)
        % load the data
        load(strcat(dataDIR,sessionName))
        
        
        jrgeco = data.jrgeco;
        chbt = data.chbo + data.chbr;    
    
        offset = runningBlocks(run).offset;

        % Initial Rest
        a = offset+startFrame1;
        b = offset+startFrame1+2*nframes;
        
        % Sustained rest
        c = offset+startFrame2;
        d = offset+startFrame2+2*nframes;

        nroi = size(jrgeco, 1);
        % Handle initial rest
        ca_acwdr = zeros(nroi, 1);
        ca_acw0 = zeros(nroi, 1);
        
        for i = 1:nroi
            [ca_acwdr(i), ca_acw0(i), ~] = acw_f(jrgeco(i,a:b), fs);
        end

        runningBlocks(run).initrest_acw.ca_acwdr = ca_acwdr;
        runningBlocks(run).initrest_acw.ca_acw0 = ca_acw0;
        runningBlocks(run).initrest_acw.ca_acw50 = ca_acw50;
        
        % now do the sustained rest
        ca_acwdr = zeros(nroi, 1);
        ca_acw0 = zeros(nroi, 1);
        
        for i = 1:nroi
            [ca_acwdr(i), ca_acw0(i), ~] = acw_f(jrgeco(i,c:d), fs);
        end
        
        runningBlocks(run).susrest_acw.ca_acwdr = ca_acwdr;
        runningBlocks(run).susrest_acw.ca_acw0 = ca_acw0;        
        
        clear jrgeco chbt
        clear rho_jRGECO rho_HbT
        clear offset
        
    end
    
    %
    save(strcat(acwDIR,'post_offset/',mousename),'runningBlocks','-v7.3')
    
end