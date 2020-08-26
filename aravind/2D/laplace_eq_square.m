%% Solve Laplace Equation on 2D plane
%% construct mesh
f = 200;
x = 0:1/f:1;
y = 0:1/f:1;

[X,Y] = meshgrid(x,y);

V = [X(:), Y(:)];
F = delaunay(V);

%% collect boundary condition points
B = unique(reshape(outline(F),[],1));
B1 = B(V(B,1) == 0); % x = 0
B2 = B(V(B,1) == 1); % x = 1
%B3 = V((V(:,1)-0.5).^2 + (V(:,2)-0.5).^2 <= 0.2^2 , :); % (x,y) in circle of radius 0.2 centered at (0.5,0.5)
B3 = find((V(:,1)-0.5).^2 + (V(:,2)-0.5).^2 <= 0.2^2);
B = [B1; B2; B3];

% define boundary condition on boundary condition points
% g = [g(BP1); g(BP2); g(BP3)]
g1 = -V(B1,2).^2;
g2 = V(B2,2).^2;
g3 = ones(length(B3), 1);
g = [g1; g2; g3];

%% solve
u = laplace_eq_2D(V, F, [B, g]);

%% plot
trisurf(F, V(:,1), V(:,2), zeros(size(V(:,1))), u, 'LineStyle','none');
%t = tsurf(F,[V u],fphong, falpha(1,0));

colormap(cbrewer('RdYlBu',40));
colorbar();

title('\nabla^2u = 0, with u = -y^2 at x = 0, u = y^2 at x = 1, u = 1 inside circle of r = 0.2 at center');
xlabel('x \in [0,1]');
ylabel('y \in [0,1]');

axis equal;
axis off;
grid off;
view(0,90);