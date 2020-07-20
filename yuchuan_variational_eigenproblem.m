%% Input data
a = 0; %left endpoint
b = 1; %right endpoint
n = 1000; %number of intervals/size of partition
h = (b-a)/n; %interval length (assume evenly-partitioned)



%% Compute Laplacian matrix L
zero_col = zeros(n,1);
G = (1/h)*([diag(ones(1,n)) zero_col] + [zero_col diag(-1*ones(1,n))]);
A = diag(h*ones(1,n)); %Area matrix (diagonal entries are the lengths of each interval)
L = (G.')*A*G; %Laplacian matrix



%% Compute eigenvalues and eigenvectors of L

[V,D] = eig(L);
eigvalues = diag(D);
min_eigvalue = eigvalues(1);
eigfunction = (1/sqrt(h))*V(:,1);




%% Plot graph of y = f(x)

x = linspace(a,b,n+1);



%% Plot first 5 eigenfunctions

for i = 1:5
    plot (x,(1/sqrt(h))*V(:,i));
    axis([0 1 -2 2])
    hold on
end

legend('eigenfunction 1','eigenfunction 2','eigenfunction 3','eigenfunction 4','eigenfunction 5');
