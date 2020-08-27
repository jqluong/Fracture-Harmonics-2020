function v = segment_metric(x, f, g)
    % function to compute a simple metric on two segment based functions.
    % computes the L2 norm of the difference function (trapezoid rule). 
    % for functions following the convention:
    %     [ y_{1,2}; y_{2,3}; ... ; y_{n,n+1}; y_{2,1}; y_{3,2}; ... ;
    %     y_{n+1,n} ]
    % takes 3 parameters: 
    %   x is function x values, f and g are functions to compare
    v = segment_norm(f-g,2,x);
end