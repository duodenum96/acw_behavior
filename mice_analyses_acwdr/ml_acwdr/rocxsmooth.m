%% Helper function
function rocx_mod = rocxsmooth(roc_x)
    rocx_mod = roc_x;
    for i = 2:length(rocx_mod)
        if rocx_mod(i) <= rocx_mod(i-1)
            rocx_mod(i) = rocx_mod(i) + abs(rocx_mod(i) - rocx_mod(i-1)) + eps;
        end
    end
end