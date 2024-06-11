function l = legend_sw(colors, location, varargin)
    %% Draw legend for swarmchart (just for convenience)
    mice = ["cm124"    "cm125"    "cm126"    "cm127"   "cm128"];
    if isempty(colors) || (nargin < 1)
        colors = [102,194,165
                  252,141,98
                  141,160,203
                  231,138,195
                  166,216,84] ./ 255;
    end
    if nargin < 2
        location = nan;
    end
    
    hold on
    L1 = scatter(nan, nan, [], colors(1, :), 'filled');
    L2 = scatter(nan, nan, [], colors(2, :),'filled');
    L3 = scatter(nan, nan, [], colors(3, :),'filled');
    L4 = scatter(nan, nan, [], colors(4, :),'filled');
    L5 = scatter(nan, nan, [], colors(5, :),'filled');
    if isnan(location)
        l = legend([L1, L2, L3, L4, L5], mice, 'box', 'off');
    else
        l = legend([L1, L2, L3, L4, L5], mice, 'box', 'off', 'Location', location);
    end
    end