function colmat = colorgenerator(colors, indx)
    %% Generate RGB colors for every category in the indx
    % Index is a vector that goes 1111222222223344444444444
    
    if isempty(colors)
        colors = [102,194,165
              252,141,98
              141,160,203
              231,138,195
              166,216,84] ./ 255;
    end
    
    ncol = size(indx, 1);
    colmat = zeros(ncol, 3);
    n_colors = length(unique(indx));
    for i = 1:n_colors
        index = indx == i;
        colmat(index, :) = repmat(colors(i, :), sum(index), 1);
    end
    end