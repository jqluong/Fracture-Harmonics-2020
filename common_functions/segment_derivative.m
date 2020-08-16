function G = segment_derivative(n, h)
    % discrete derivative for the segment representation of a function.
    % for functions following the convention:
    %     [ y_{1,2}; y_{2,3}; ... ; y_{n,n+1}; y_{2,1}; y_{3,2}; ... ;
    %     y_{n+1,n} ]
    % takes 2 parameters: n is number of segments, h segment width.
    % h defaults to 1.
    % output is matrix G such that G*y = [y_1'; y_2'; ... ; y_n']
    
    if nargin > 1 % check if h is supplied, else default to 1
        w = h;
    else
        w = 1;
    end
    
    G = 1/w * kron([-1 1], speye(n));
end