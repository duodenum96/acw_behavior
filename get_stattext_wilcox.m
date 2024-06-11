function stattext = get_stattext_wilcox(z, p, r, orientation)

    if (p > 0.001)
        if string(orientation) == "vertical"
            pstring = sprintf('p = %.3f', p);
        elseif string(orientation) == "horizontal"
            pstring = sprintf('p = %.3f, ', p);
        end
    elseif p <= 0.001
        if string(orientation) == "vertical"
            pstring = 'p < 0.001';
        elseif string(orientation) == "horizontal"
            pstring = 'p < 0.001, ';
        end
    end
    
    if string(orientation) == "vertical"
        stattext = {sprintf('z = %.2f', abs(z)), ...
            pstring, ...
            sprintf('r = %.2f', r)};
    elseif string(orientation) == "horizontal"
        stattext = [sprintf('z = %.2f, ', abs(z)), ...
            pstring, ...
            sprintf('r = %.2f', r)];
    end
end