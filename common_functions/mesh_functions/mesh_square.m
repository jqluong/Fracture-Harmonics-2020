
function [V,F] = mesh_square(s)
% generates square mesh
x = 0:1:s;
y = x';
[x,y] = meshgrid(x,y);
x = x(:);
y = y(:);
F = delaunay(x,y);
V = [x y];
end