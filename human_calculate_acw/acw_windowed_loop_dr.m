function [acwdrmat] = acw_windowed_loop_dr(ts, fs, windowsize)
    %% Loop over channels
    % ts is a matrix (channels x time)
    % fs is sampling frequency in Hz
    % windowsize is length of window in seconds
    
    nchan = size(ts, 1);
    acwdrcell = cell(nchan, 1);
    for i = 1:nchan
        acwdrcell{i} = acw_windowed_dr(ts(i, :), fs, windowsize);
    end
    
    acwdrmat = cat(1, acwdrcell{:});
    end