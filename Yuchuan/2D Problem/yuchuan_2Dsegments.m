[x, y] = meshgrid(1:10, 1:10);
points = [x(:), y(:)];
V = [points zeros(size(points,1),1)];
F = delaunay(points);

u = ones(3*size(F,1),1);
u(1:2) = 8;

D = discontinuity_matrix(V,F);
Du = sparse(D*u);