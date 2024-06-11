function runningBlocks = running_acws(mousename,runningBlocks,acwDIR,dataDIR,delay, fs)
    % This function computes the correlation map over the running window for each
    % running block. The running window is defined as 10 s middle of locomotion.
    % For each mouse, it returns a struct of running blocks with two new fields: 
    % correlations_neural and correlations_HbT, which contain correlation maps
    % over the running winow for neural and hemo data.
    % 
    % Author: Somayeh "Bahar" Shahsavarani
    % email: bahar@huskers.unl.edu
    % 
    % Modified by Yasir for ACW values
         
    % remove running bouts with less than 20 s duration
    duration = [runningBlocks.duration];
    durationIDX = duration < 400;
    runningBlocks(durationIDX) = [];
    
    % for FC maps of 10 seconds, set it to 100 for half a window
    nframes = 100;
    
    dt = 1 / fs;
    for run = 1:length(runningBlocks)
        sessionName = runningBlocks(run).session;
        
        fprintf('run %d / %d \n', run, length(runningBlocks))
        % load the data
        load(strcat(dataDIR,sessionName))
        
        jrgeco = data.jrgeco;
        chbt = data.chbo + data.chbr;
            
        duration = runningBlocks(run).duration;
        onset = runningBlocks(run).onset;
        
        a = onset+floor(duration/2)-nframes;
        b = onset+floor(duration/2)+nframes;
        
        nroi = size(jrgeco, 1);
        ca_acwdr = zeros(nroi, 1);
        ca_acw0 = zeros(nroi, 1);
        
        for i = 1:nroi
            [ca_acwdr(i), ca_acw0(i), ca_acw50(i)] = acw_f(jrgeco(i,a:b), fs);
        end
        
        runningBlocks(run).running_acw.ca_acwdr = ca_acwdr;
        runningBlocks(run).running_acw.ca_acw0 = ca_acw0;
        
        
        clear jrgeco chbt
        clear rho_jRGECO rho_HbT
        clear onset       
    end
    save(strcat(acwDIR,'running/',mousename),'runningBlocks','-v7.3')
    
end
