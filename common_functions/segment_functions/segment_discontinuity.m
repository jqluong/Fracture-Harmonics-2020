function D = segment_discontinuity(n)
    % matrix to compute discontinuity between segments
    % for functions following the convention:
    %     [ y_{1,2}; y_{2,3}; ... ; y_{n,n+1}; y_{2,1}; y_{3,2}; ... ;
    %     y_{n+1,n} ]
    % takes 1 parameters: n is number of segments
    % output is a matrix that produces a vector with all of the 
    % discontinuities as elements:
    %     [ y_{2,3}-y_{2,1}; y_{3,4}-y_{3,2}; ... ; y_{n,n+1}-y_{n,n-1} ]
    
    D = [sparse(1:n-1,2:n,1,n-1,n),-speye(n-1,n)];
end