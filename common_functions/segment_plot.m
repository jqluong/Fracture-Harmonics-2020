function segment_plot(x, y, c)
    % this function will plot a discrete function defined by x and y, in
    % color c.
    % for functions following the convention (same for x):
    %     [ y_{1,2}; y_{2,3}; ... ; y_{n,n+1}; y_{2,1}; y_{3,2}; ... ;
    %     y_{n+1,n} ]
    % takes three parameters: x is x values, y is y values, c is color
    % pass c in as: 'b', or 'r'. defaults to 'b'
    
    if nargin == 2 % check if h is supplied, else default to 1
        c = 'b';
    end
    
    n = length(x)/2;
    
    x1 = x(1:n)'; % line segment position 1
    x2 = x(n+1:end)'; % line segment position 2
    
    y1 = y(1:n)';
    y2 = y(n+1:end)';
    
    plot([x1;x2],[y1;y2], 'Color', c);
end