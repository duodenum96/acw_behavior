function acwdrs = acw_windowed_dr(ts, fs, windowsize)
    %% Calculate ACWs in sliding window fashion with no overlap
    % ts is 1D vector
    % Windowsize is in unit of seconds
    % fs is sampling frequency in Hz
    
    window_samples = windowsize * fs;
    datalength = length(ts);
    
    acwdrs = [];
    i = 0;
    while true
        indx = ((i*window_samples)+1):((i+1)*window_samples);
        if indx(end) > datalength
            break
        end
    
        acwdr = acw_dr(ts(indx), fs);
        acwdrs = [acwdrs, acwdr];
        i = i + 1;
    end
    end
    