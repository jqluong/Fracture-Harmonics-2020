
function [V,F] = mesh_square(s)
% generates square mesh
[x, y] = meshgrid(1:s, 1:s);
pts = [x(:), y(:)];
V = [pts zeros(size(pts,1),1)];
F = delaunay(pts);
end