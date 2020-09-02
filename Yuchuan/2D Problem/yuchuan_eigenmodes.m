[x, y] = meshgrid(1:10, 1:10);
points = [x(:), y(:)];
V = [points zeros(size(points,1),1)];
F = delaunay(points);

u = zeros(3*size(F,1),1);
u(1)=1;
u(2) = 1;
u(11) = 1;



G = face_grad(V,F);
Gu = G*u;

face_plot(V,F,u);