%jack's testing script, it might get moved to my folder later
%works with normal old face functions
x = 0:.1:1;
y = x';
[x,y] = meshgrid(x,y);
x = x(:);
y = y(:);
f = x + y;
F = delaunay(x,y);
V = [x y];
fface = vertexToFace(F, f);
face_plotting(V, F, fface)
%works with discontinous functions
clear all
x = [0 0 1 1]';
y = [0 1 0 1]';
F = delaunay(x,y);
V = [x y];
f = [0 0 0 1 1 1]';
face_plotting(V,F,f)