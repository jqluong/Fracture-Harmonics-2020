%% problem
a = 0; % a, b are endpoints
b = 1;
fa = 0; % function value at endpoints
fb = 1;
n = 101; % number of points
h = (b-a)/(n-1); % step size

%% solve
% first, generate appropriate laplacian matrix
L = discrete_laplacian(n,h);
f = zeros(n,1); % solution

f(1) = fa;
f(n) = fb;

r = -L*f;
f(2:n-1) = L(2:n-1 , 2:n-1) \ r(2:n-1);

plot(0:h:1, f);