function runningBlocks = onset_acws(mousename,runningBlocks,acwDIR,dataDIR,delay,fs)
    % This function computes the correlation map over the onset window for each
    % running block. The onset window overlaps locomotion and rest with 5 s
    % before locomotion and 5 s beginning of locomotion.
    % For each mouse, it returns a struct of running blocks with two new fields: 
    % correlations_neural and correlations_HbT, which contain correlation maps
    % over the onset winow for neural and hemo data.
    % 
    % Author: Somayeh "Bahar" Shahsavarani
    % email: bahar@huskers.unl.edu
    % Modified by Yasir for ACW values
    
    % remove running bouts with less than 10 s duration
    duration = [runningBlocks.duration];
    durationIDX = duration < 200;
    runningBlocks(durationIDX) = [];
    
    % remove running bouts with pre-running rest time less than 60 s
    time2previousrun = [runningBlocks.time2previousrun];
    time2previousrunIDX = time2previousrun < 1200;
    runningBlocks(time2previousrunIDX) = [];
    
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
            
        onset = runningBlocks(run).onset;
        
        a = onset-nframes;
        b = onset+nframes;

        nroi = size(jrgeco, 1);
        ca_acwdr = zeros(nroi, 1);
        ca_acw0 = zeros(nroi, 1);
        
        for i = 1:nroi
            [ca_acwdr(i), ca_acw0(i), ~] = acw_f(jrgeco(i,a:b), fs);
        end
        
        runningBlocks(run).onset_acw.ca_acwdr = ca_acwdr;
        runningBlocks(run).onset_acw.ca_acw0 = ca_acw0;
        
        
        clear jrgeco chbt
        clear rho_jRGECO rho_HbT
        clear onset
        
    end
    
    %
    save(strcat(acwDIR,'onset/',mousename),'runningBlocks','-v7.3')
    
end