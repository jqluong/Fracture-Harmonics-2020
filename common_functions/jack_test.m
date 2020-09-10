%jack's testing script, it might get moved to my folder later
%works with normal old face functions
mkdir('vertex_functions')
addpath('vertex_functions')
x = 0:.1:1;
y = x';
[x,y] = meshgrid(x,y);
x = x(:);
y = y(:);
F = delaunay(x,y);
V = [x y];
u = laplace_eq_2D_seg_quadprog(V,F);
face_plotting(V,F,u)