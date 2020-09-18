%% first, construct mesh
f = 40;
x = -2:1/f:2;
y = -2:1/f:2;

[X,Y] = meshgrid(x,y);

V = [X(:), Y(:)];
F = delaunay(V);

%% construct outline
%n = 1024; % boundary spacing: 2pi/n
%th = [0:2*pi()/n:2*pi()]';
%x = cos(th);
%y = sin(th);

%% constuct mesh
%[V, F] = triangle([x,y], 'Quality', 30, 'MaxArea', 0.01);

%% construct laplacian
%G = face_grad(V, F);
G = grad(V,F);
L = G'*G;

%% find eigenvalues
[E,D] = eigs(L);

%% plot
%face_plot(V, F, E(:,1));
%face_plot(V, F, E(:,2));
%face_plot(V, F, E(:,3));
trisurf(F, V(:,1), V(:,2), E(:,2), 'LineStyle','none');
colormap(cbrewer('RdYlBu',40));
colorbar();