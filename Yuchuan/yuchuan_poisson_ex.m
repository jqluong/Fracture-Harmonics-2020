%% Example of solving Possion equation with Dirichlet boundary condition
%requires poisson_2D.m

clear


%% Generate mesh
[x, y] = meshgrid(1:10, 1:10);
points = [x(:), y(:)];
V = [points zeros(size(points,1),1)];
F = delaunay(points);


%% Input
n = size(V,1);
f = sparse(n,1); %source term in Poisson equation

%% Create vector for boundary values u_0
O = outline(F);     %boundary outline (edges)
bv = unique(O(:));  %boundary vertices
v = zeros(n,1);

for i = bv
    v(i) = 100;
end 

u_0 = v;


%% Solve for u
u_int = poisson_2D(V,F,f,u_0);
u = u_int + u_0;


