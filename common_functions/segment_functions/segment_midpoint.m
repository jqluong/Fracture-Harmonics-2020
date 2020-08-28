function M = segment_midpoint(n)
    % function that returns matrix to find the midpoint of each segment in
    % a segment based discrete function
    % for functions following the convention:
    %     [ y_{1,2}; y_{2,3}; ... ; y_{n,n+1}; y_{2,1}; y_{3,2}; ... ;
    %     y_{n+1,n} ]
    % takes 1 parameters: n is number of segments
    % output is matrix M such that M*y = 1/2 * [ y_{1,2}+y_{2,1}; ... ]
    
    M = 1/2 * kron([1 1], speye(n));
end