%% Solve Laplace Equation on 2D plane

% first, construct mesh
f = 200;
x = 0:1/f:1;
y = 0:1/f:1;

[X,Y] = meshgrid(x,y);

V = [X(:), Y(:)];
F = delaunay(V);

% collect boundary condition points, nice to use inequality
BP1 = V(V(:,1) == 0,:); % x = 0
BP2 = V(V(:,1) == 1,:); % x = 1
BP3 = V((V(:,1)-0.5).^2 + (V(:,2)-0.5).^2 <= 0.2^2 , :); % (x,y) in circle of radius 0.2 centered at (0.5,0.5)
B = [BP1; BP2; BP3];

% define boundary condition on boundary condition points
% g = [g(BP1); g(BP2); g(BP3)]
g = [-BP1(:,2).^2; BP2(:,2).^2; ones(length(BP3), 1)];
u = laplace_eq_2D(V, F, [B, g]);

trisurf(F, V(:,1), V(:,2), zeros(size(V(:,1))), u, 'LineStyle','none');
%trisurf(F, V(:,1), V(:,2), u,'LineStyle','none');
colorbar();
title('\nabla^2u = 0, with u = -y^2 at x = 0, u = y^2 at x = 1, u = 1 inside circle of r=0.2 at center');
xlabel('x \in [0,1]');
ylabel('y \in [0,1]');