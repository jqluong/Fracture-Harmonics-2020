
function [V,F] = mesh_square(s)
% generates square mesh
% s: number of partitions on interval from 0 to 1
x = [0:s:1];
y = x';
[x,y] = meshgrid(x,y);
x = x(:);
y = y(:);

F = delaunay(x,y);
V = [x y];
end