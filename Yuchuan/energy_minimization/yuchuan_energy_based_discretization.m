%% Input data
a = 0; %left endpoint
b = 1; %right endpoint
n = 100; %number of intervals/size of partition
h = 1/n; %interval length (assume evenly-partitioned)
f_a = 0; %left boundary condition f(a) = 0
f_b = 1; %right boundary condition f(b) = 1



%% Compute Laplacian matrix L
zero_col = zeros(n,1);
G = (1/h)*([diag(ones(1,n)) zero_col] + [zero_col diag(-1*ones(1,n))]);
A = diag(h*ones(1,n)); %Area matrix (diagonal entries are the lengths of each interval)
L = (G.')*A*G; %Laplacian matrix



%% Solve for energy-minimizing f (with given boundary conditions)

f_endpoints = zeros(n+1,1); 
f_endpoints(1) = f_a;
f_endpoints(n+1) = f_b;
c = -L*f_endpoints;
f_interior = L(2:n , 2:n)\c(2:n);
f = [f_a;f_interior;f_b];




%% Plot graph of y = f(x)

x = linspace(a,b,n+1);
plot(x,f);






