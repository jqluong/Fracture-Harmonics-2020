function L = discrete_derivative(n,h)
    L = zeros(n,n);
    if n < 3
        error("too few points");
    else
        L(1,1) = -3;
        L(1,2) = 4;
        L(1,3) = -1;
        L(n,n-2) = 1;
        L(n,n-1) = -4;
        L(n,n) = 3;

    for i = 2:n-1
        L(i,i-1) = -1;
        L(i,i+1) = 1;
    end
    L = 1/(2*h)*L;
end