%% Example of solving Possion equation with Dirichlet boundary condition
%requires poisson_2D.m

clear


%% Generate mesh
[x, y] = meshgrid(1:10, 1:10);
points = [x(:), y(:)];
V = [points zeros(size(points,1),1)];
F = delaunay(points);
z = zeros(size(x));



%% Input
n = size(V,1);
f = sparse(n,1); %source term in Poisson equation


E = edges(F);


%% 1. Boundary condition problem

% % Create vector for boundary values u_0
% O = outline(F);     %boundary outline (edges)
% bv = unique(O(:));  %boundary vertices
% v = zeros(n,1);
% 
% % for i = bv               %%this is for the boundary condition of u(x,y) = x + y along the boundary
% %     v(i) = V(i,1)+V(i,2);
% % end 
% 
% u_0 = v;
% 
% % Solve for u
% u_int = poisson_2D(V,F,f,u_0,bv);
% v_int = setdiff((1:size(V,1)),bv).'; %interior vertices
% s = sparse(size(V,1),1);
% s(v_int) = u_int;
% u = s + u_0;





%% 2. Unit norm constraint problem

u = poisson_2D(V,F,f);







%% Plot u


X = reshape(V(:,1),10,10)';
Y = reshape(V(:,2),10,10)';

%for boundary condition problem
%Z = reshape(u,10,10)';   
%surf(X,Y,Z);

%for unit norm constraint, plots first 6 eigenmodes
figure
subplot(2,3,1)
surf(X,Y,reshape(u(:,1),10,10)')
subplot(2,3,2)
surf(X,Y,reshape(u(:,2),10,10)')
subplot(2,3,3)
surf(X,Y,reshape(u(:,3),10,10)')
subplot(2,3,4)
surf(X,Y,reshape(u(:,4),10,10)')
subplot(2,3,5)
surf(X,Y,reshape(u(:,5),10,10)')
subplot(2,3,6)
surf(X,Y,reshape(u(:,6),10,10)')
