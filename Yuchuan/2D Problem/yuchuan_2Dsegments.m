[x, y] = meshgrid(1:10, 1:10);
points = [x(:), y(:)];
V = [points zeros(size(points,1),1)];
F = delaunay(points);

u = zeros(3*size(F,1),1);
u(1:3)=1;

D = discontinuity_matrix(F);

Du = sparse(D*u);