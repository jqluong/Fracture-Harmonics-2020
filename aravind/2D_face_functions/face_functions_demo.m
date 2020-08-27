% program to test face defined function programs
%% first, construct mesh
f = 10;
x = -2:1/f:2;
y = -2:1/f:2;

[X,Y] = meshgrid(x,y);

V = [X(:), Y(:)];
F = delaunay(V);

%% construct face function u = x^2 - y^3
u = zeros(3*size(F,1), 1);
for i = 1:length(F)
    for j = [2 1 0]
        u(3*i - j) = V(F(i, 3-j), 1)^2 - V(F(i, 3-j), 2)^3; % x^2 - y^3
    end
end

%% find gradient of u
G = face_grad(V, F); % discrete gradient operator
du = G*u;
% du is [dx_1 dy_1 ... dx_n, dy_n1]'
du = [du(1:2:end), du(2:2:end)];
% reshape(du, size(V,2), [])'

%% find centers of each face
C = face_centers(V, F);
%scatter(C(:,1), C(:,2))

%% plot gradient
quiver(C(:,1), C(:,2), du(:,1), du(:,2))
axis equal

%% compute Dirichlet Energy
E = face_dirichlet_energy(V, F, u, G)

%% plot u
face_plot(V, F, u);