clear all 
opengl software
p = 10;
xaxis = 0:1/p:1;
yaxis = 0:1/p:1;
[xaxis, yaxis] = meshgrid(xaxis, yaxis);
xaxis = xaxis(:);
yaxis = yaxis(:);
[numberOfVertices, ~] = size(xaxis);
T = delaunay(xaxis,yaxis);
G = grad([xaxis yaxis], T);
cellArea = doublearea([xaxis yaxis], T)/2;
[numberOfFaces, ~] = size(cellArea);
cellAreaMatrix = [diag(cellArea) sparse(numberOfFaces, numberOfFaces); sparse(numberOfFaces, numberOfFaces) diag(cellArea)];
constraintMatrix = [speye(p) sparse(p, numberOfVertices-p); sparse(numberOfVertices-2*p, numberOfVertices);sparse(p, numberOfVertices-p) speye(p)];
constraintVector = [zeros(numberOfVertices - p, 1); ones(p,1)];
options = optimoptions(@quadprog, 'Algorithm', 'interior-point-convex', 'OptimalityTolerance', 1e-12);
quadraticTerm = G'*cellAreaMatrix*G;
f = quadprog(quadraticTerm/2, zeros(numberOfVertices,1) , zeros(numberOfVertices,numberOfVertices), zeros(numberOfVertices,1),constraintMatrix, constraintVector);
trisurf(T, xaxis, yaxis, f)