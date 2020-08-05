%% problem
a = 0; % endpoints
b = 1;
n = 100;
h = (b-a)/(n-1);

% generate laplacian
L = discrete_laplacian(n, h);

% laplacian eigenvalues/eigenfunctions
[F,V] = eigs(L);
F = 1/sqrt(h) * F;

% plot eigenfunctions
hold on;
for i = 1:5
   plot(a:h:b, F(:,i)) 
end