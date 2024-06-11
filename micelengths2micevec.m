function micevec = micelengths2micevec(micelenghts)
    %% Helper function to turn a vector of lengths to a vector of 1s, 2s etc for each length
    
    micevec = [];
    for i = 1:length(micelenghts)
        micevec = [micevec; ones(micelenghts(i), 1)*i];
    end
    end