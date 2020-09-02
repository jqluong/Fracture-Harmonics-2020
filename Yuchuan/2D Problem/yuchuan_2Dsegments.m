[x, y] = meshgrid(1:10, 1:10);
points = [x(:), y(:)];
V = [points zeros(size(points,1),1)];
F = delaunay(points);

u = zeros(3*size(F,1),1);
u(1:2)=1;
u(15) = 1;




face_plotting(V,F,u);


D = discontinuity_matrix(F);
Du = sparse(D*u);

G = face_grad(V,F);
Gu = sparse(G*u);

M = face_area_matrix(V,F);