%Solving Laplace's equation on unit square with boundary conditions f(0,y)
%= 0 and f(1,y) = 1.  Did so by minimization of Dirichlet Energy via
%quadprog on a triangle mesh.
clear all 
opengl software
p = input('Please give number of nodes: ') - 1;
%Creates the axes for triangular mesh
xaxis = 0:1/p:1;
yaxis = 0:1/p:1;
[xaxis, yaxis] = meshgrid(xaxis, yaxis);
xaxis = xaxis(:);
yaxis = yaxis(:);
[numberOfVertices, ~] = size(xaxis);
%Creates face matrix and gradient based on triangular mesh
T = delaunay(xaxis,yaxis);
G = grad([xaxis yaxis], T);
%Finds area of each face and creates the area matrix used in discretization
%of Dirichlet Energy
cellArea = doublearea([xaxis yaxis], T)/2;
[numberOfFaces, ~] = size(cellArea);
cellAreaMatrix = [spdiags(cellArea, [0], numberOfFaces, numberOfFaces) sparse(numberOfFaces, numberOfFaces); sparse(numberOfFaces, numberOfFaces) spdiags(cellArea, [0], numberOfFaces, numberOfFaces)];
quadraticTerm = G'*cellAreaMatrix*G;
%Converts boundary conditions
constraintMatrix = [speye(p) sparse(p, numberOfVertices-p); sparse(numberOfVertices-2*p, numberOfVertices);sparse(p, numberOfVertices-p) speye(p)];
constraintVector = [sparse(numberOfVertices - p, 1); ones(p,1)];
%Solve minimization problem and plots the results
options = optimoptions(@quadprog, 'Algorithm', 'interior-point-convex', 'OptimalityTolerance', 1e-12);
f = quadprog(quadraticTerm/2, zeros(numberOfVertices,1) , [], [],constraintMatrix, constraintVector);
trisurf(T, xaxis, yaxis, f)
title('Solve Laplace Equation with BC f(0,y) = 0, f(1,y) = 1')