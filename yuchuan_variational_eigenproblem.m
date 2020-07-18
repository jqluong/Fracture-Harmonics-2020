%% Input data
a = 0; %left endpoint
b = 1; %right endpoint
n = 10; %number of intervals/size of partition
h = 1/n; %interval length (assume evenly-partitioned)
f_a = 0; %left boundary condition f(a) = 0
f_b = 1; %right boundary condition f(b) = 1



%% Compute Laplacian matrix L
zero_col = zeros(n,1);
G = (1/h)*([diag(ones(1,n)) zero_col] + [zero_col diag(-1*ones(1,n))]);
A = diag(h*ones(1,n)); %Area matrix (diagonal entries are the lengths of each interval)
L = (G.')*A*G; %Laplacian matrix



%% Compute eigenvalues and eigenvectors of L

[V,D] = eig(L);
eigvalues = diag(D);
min_eigvalue = eigvalues(1);
eigfunction = V(:,1);



%% Plot graph of y = f(x)

x = linspace(a,b,n+1);
plot(x,eigfunction);


