function segment_plot_legend(l, c)
    % function to plot a color coded legend.
    % in matlab, legend is tied to the actual plot, so need to work around
    % to have legends for these segment functions
    % takes 2 inputs: l is string array for legend titles, c is array of
    % corresponding colors
    % note curly brace needed for string array, square for color
    %   eg. l = {'func1', 'func2', ... }, c = ['r', 'b', ...]
    
    % plot point at NaN - won't be displayed, so can tie the colors and
    % legend to that
    h = zeros(size(l));
    for i=1:length(h)
        h(i) = plot(NaN,NaN,c(i)); % plot point at NaN - won't be displayed
    end
    legend(h, l);
end