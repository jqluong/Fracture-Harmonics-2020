function L = discrete_laplacian(n,h)
    L = zeros(n,n);
    if n < 3
        error("too few points");
    elseif n == 3
        L(1,1) = 1;
        L(1,2) = -2;
        L(1,3) = 1;
        L(n,n-2) = 1;
        L(n,n-1) = -2;
        L(n,n) = 1;
    else
        L(1,1) = 2;
        L(1,2) = -5;
        L(1,3) = 4;
        L(1,4) = -1;
        L(n,n-3) = -1;
        L(n,n-2) = 4;
        L(n,n-1) = -5;
        L(n,n) = 2;
    end
    for i = 2:n-1
        L(i,i-1) = 1;
        L(i,i) = -2;
        L(i,i+1) = 1;
    end
    L = 1/h^2*L;
end