%% solve Laplace equation on annulus
%% construct outline
% get points on outer ring of annulus
ro = 1;
th = [0:2*pi()/1024:2*pi()]';
x = ro*cos(th);
y = ro*sin(th);
V = [x, y];
P = [1:length(th)]';
E = [P, P+1]; % construct edges connecting outer ring vertices
E(end,2) = 1;
n1 = length(E);

% get points on inner ring of annulus
ri = ro/4;
th = [0:2*pi()/512:2*pi()]';
x = ri*cos(th);
y = ri*sin(th);
V = [V; x, y];
P = [1:length(th)]' + n1;
E = [E; P, P+1]; % construct edges connecting inner ring vertices
E(end,2) = 1 + n1;

%% construct mesh
% generate mesh, specifing [0,0] as located in a hole
[V, F] = triangle(V, E, [0,0],'Quality', 30, 'MaxArea', 0.001);
%triplot(F,V(:,1),V(:,2));

%% process boundary conditions
% Boundary conditions
% on outer boundary, g = sin(4*theta)
% on inner boundary, g = -1, y < 0
%                    g = 1, y > 0

% collect points on boundary
B = unique(reshape(outline(F),[],1));

% determine inner boundary or outer boundary
B1 = B(abs(normrow(V(B,:)) - ro) <= 0.1); % outer
B2 = B(abs(normrow(V(B,:)) - ri) <= 0.1); % inner
B2p = B2(V(B2, 2) >= 0); % inner; y >= 0
B2m = B2(V(B2, 2) < 0); % inner; y < 0
B = [B1; B2p; B2m];

g1 = sin(4*(atan(V(B1,1) ./ V(B1,2))+pi()/8));
g2p = ones(length(B2p),1);
g2m = -ones(length(B2m),1);
g = [g1; g2p; g2m];

%% solve Laplace eq
u = laplace_eq_2D(V, F, [B, g]);

%% plot
%trisurf(F, V(:,1), V(:,2), zeros(size(V(:,1))), u, 'LineStyle','none');
t = tsurf(F,[V u],fphong, falpha(1,0));

colormap(cbrewer('RdYlBu',40));
colorbar();

title('\nabla^2u = 0, with u = cos(4\theta) on outer boundary, u = \pm 1 on inner boundary');

axis equal
grid off;
axis off;
camlight;
add_isolines(t);