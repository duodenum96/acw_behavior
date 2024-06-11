function runningBlocks = offset_acws(mousename,runningBlocks,acwDIR,dataDIR,delay,fs)
    % This function computes the correlation map over the offset window for each
    % running block. The offset window overlaps locomotion and rest with 5 s
    % end of locomotion and 5 s beginning of rest.
    % For each mouse, it returns a struct of running blocks with two new fields: 
    % correlations_neural and correlations_HbT, which contain correlation maps
    % over the offset winow for neural and hemo data.
    % 
    % Author: Somayeh "Bahar" Shahsavarani
    % email: bahar@huskers.unl.edu
    % Modified by Yasir
    
    % remove running bouts with less than 20 s duration
    duration = [runningBlocks.duration];
    durationIDX = duration < 400;
    runningBlocks(durationIDX) = [];
    
    % remove running bouts with post-running rest time less than 10 s
    time2nextrun = [runningBlocks.time2nextrun];
    time2nextrunIDX = time2nextrun < 200;
    runningBlocks(time2nextrunIDX) = [];
    
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
        
        a = offset-nframes;
        b = offset+nframes;
        
        nroi = size(jrgeco, 1);
        ca_acwdr = zeros(nroi, 1);
        ca_acw0 = zeros(nroi, 1);
        
        for i = 1:nroi
            [ca_acwdr(i), ca_acw0(i), ~] = acw_f(jrgeco(i,a:b), fs);
        end
        
        runningBlocks(run).offset_acw.ca_acwdr = ca_acwdr;
        runningBlocks(run).offset_acw.ca_acw0 = ca_acw0;
        
        
        clear jrgeco chbt
        clear rho_jRGECO rho_HbT
        clear offset
        
    end
    
    %
    save(strcat(acwDIR,'offset/',mousename),'runningBlocks','-v7.3')
    
end