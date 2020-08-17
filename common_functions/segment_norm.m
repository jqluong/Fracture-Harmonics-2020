function v = segment_norm(y, p, x)
    % function to compute the Lp norm of input segment based function.
    % uses the midpoint of the segment as the segment value; trapezoid rule. 
    % for functions following the convention:
    %     [ y_{1,2}; y_{2,3}; ... ; y_{n,n+1}; y_{2,1}; y_{3,2}; ... ;
    %     y_{n+1,n} ]
    % takes 3 parameters: 
    %   y is function, p is norm number, x is
    %   optional function x values (defaults to even spacing on [0,1])
    
    n = length(y)/2;
    
    if p == 0
        v = n;
        return;
    end
    
    if nargin == 2 % check if x is supplied, else use even segments on [0,1]
        w = 1/n * ones(n,1);
    else
        w = segment_derivative(n)*x; % width of each segment
    end
    
    ym = abs(segment_midpoint(n)*y); % midpoints of each segment
    ym = ym.^p; % raise to pth power
    v = sum(w.*ym)^(1/p);
end
