opengl software
clear all
%Dirichlet Energy
n = input('Number of Points: ');
x = 0:1/n:1;
y = x';
f = x.^2 + y.^2;
[dx, dy] = gradient(f,1/n);
gradient = dx.^2 + dy.^2;
dirichletEnergy = (1/n)^2*sum(gradient, 'all')
%Triangular mesh
xaxis = 0:1/n:1;
yaxis = xaxis;
[xaxis, yaxis] = meshgrid(xaxis, yaxis);
z = xaxis + yaxis;
T = delaunay(xaxis,yaxis);
trisurf(T, xaxis, yaxis, f)
%Dirichlet Energy on triangular mesh
xtriangle = xaxis(:);
ytriangle = yaxis(:);
G = grad([xtriangle ytriangle], T);
gradientFunction = G*f(:);
gradientFunction = gradientFunction .* gradientFunction;
dirichletEnergy2 = (1/(2*n^2))*sum(gradientFunction, 'all')