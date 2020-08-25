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

% generate mesh, specifing [0,0] as located in a hole
[V, F] = triangle(V, E, [0,0],'Quality', 30, 'MaxArea', 0.001);
%triplot(F,V(:,1),V(:,2));

% Boundary conditions
% on outer boundary, g = sin(4*theta)
% on inner boundary, g = -1, y < 0
%                    g = 1, y > 0

BP1 = V(V(:,1).^2 + V(:,2).^2 == ro^2, :);
BP2 = V(V(:,1).^2 + V(:,2).^2 == ri^2 & V(:,2) < 0, :);
BP3 = V(V(:,1).^2 + V(:,2).^2 == ri^2 & V(:,2) > 0, :);
BP = [BP1; BP2; BP3];

g1 = sin(4*(atan(BP1(:,1)./BP1(:,2))+pi()/8));
g2 = -ones(length(BP2),1);
g3 = ones(length(BP3),1);
g = [g1; g2; g3];

% solve
u = laplace_eq_2D(V, F, [BP, g]);

% plot
trisurf(F, V(:,1), V(:,2), zeros(size(V(:,1))), u, 'LineStyle','none');
colorbar();
title('\nabla^2u = 0, with u = cos(4\theta) on outer boundary, u = \pm 1 on inner boundary');
xlabel('x');
ylabel('y');