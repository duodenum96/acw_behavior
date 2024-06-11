function [acw_dr, acw_0, acw_50, acfunc, lags] = acw_f(ts, fs, trim, trimtp)
    % Translation of acw_dr_f from src.py to MATLAB
    % Parameters
    % ---------
    % ts : 1D Numpy vector (Has to be in the shape (n, ). Otherwise won't work)
    %     Time series.
    % fs : double
    %     Sampling rate (in Hz).
    % fast : Logical
    %     Use FFT to calculate ACF. The default is True.
    % trim : Logical
    %     If True, trim what comes after the ACF reaches 0.
    % trimtp : Integer
    %     Number of time points to count after ACF reaches 0 before trimming.
    %     For example, if this is 5, the ACF will be trimmed after 5 lags from ACW-0
    %     point.

    % Returns
    % -------
    % acw_dr : Autocorrelation window calculated as the decay rate of ACF (in seconds).
    % acw_0  : Autocorrelation window calculated as where ACF crosses 0.
    % acfunc : Autocorrelation function (for troubleshooting / plotting etc.)
    % lags : x-axis of ACF, for plotting purposes

    if nargin < 4
        trim = false;
        trimtp = 5;
    elseif nargin < 3
        trim = false;
    elseif nargin < 2
        error('Have to specify at least ts and fs')
    end
        
    [acfunc, lags] = autocorr(ts, length(ts)-1);

    [~, acw_0] = max(acfunc <= 0);
    [~, acw_50] = max(acfunc <= 0.5);

    lags = lags / fs;

    % Decoration: Get rid of the portion after ACW-0.
    if trim
        acfunc = acfunc(1 : (acw_0 + trimtp));
        lags = lags(1 : (acw_0 + trimtp));
    end

    % Initial guess (starting point) in nonlinear optimization.
    % Note that this is a bit arbitrary. What matters is that this has to be
    % same in all subjects
    init_guess = 0.5;
    expdecayfun = @(dr, x) exp(-(1/dr)*x);
    acw_dr = lsqcurvefit(expdecayfun, init_guess, lags, acfunc, 0, inf);
    
    acw_0 = acw_0 / fs;
    acw_50 = acw_50 / fs;
    lags = lags / fs;
end