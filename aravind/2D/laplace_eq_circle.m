%% solve Laplace equation on disk
%% construct outline
n = 1024; % boundary spacing: 2pi/n
th = [0:2*pi()/n:2*pi()]';
x = cos(th);
y = sin(th);

%% constuct mesh
[V, F] = triangle([x,y], 'Quality', 30, 'MaxArea', 0.001);

% Boundary conditions: on circle. if y < 0, u = x^3, if y > 0, u = y^2
% here, collect vertices that have boundary conditions enforced 

%% process boundary conditions
% first, collect boundary points
B = unique(reshape(outline(F),[],1));
BP1 = B(V(B,2) <= 0);
BP2 = B(V(B,2) > 0);
B = [BP1; BP2];

g1 = V(BP1, 1).^3;
g2 = V(BP2, 2).^2;
g = [g1; g2];


%% solve Laplace eq
u = laplace_eq_2D(V, F, [B, g]);

%% plot
%trisurf(F, V(:,1), V(:,2), zeros(size(V(:,1))), u, 'LineStyle','none');
t = tsurf(F,[V u],fphong, falpha(1,0));
colormap(cbrewer('RdYlBu',40));
colorbar();

title('\nabla^2u = 0, with u = x^3 at y < 0 and u = y^2 at y > 0 on unit circle ');
%xlabel('x \in [-1,1]');
%ylabel('y \in [-1,1]');

axis equal
grid off;
axis off;
camlight;
add_isolines(t);

%view(0,90);