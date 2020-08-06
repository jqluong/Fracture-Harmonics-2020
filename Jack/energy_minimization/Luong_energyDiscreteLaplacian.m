clear all
opengl software
%Ask user to input a vector of stepsizes and boundary conditions
stepsize = input("Please input a horizontal vector of stepsizes: ");
boundaryConditions0 = input("Please input the first boundary condition: ");
boundaryConditions1 = input("Please input the second boundary condition: ");

%Convert to column vectors, create a dummy function variable
stepsize = stepsize'; 
[stepnumber,~] = size(stepsize);
f = zeros(stepnumber+1, 1);
f(1) = boundaryConditions0;
f(length(f)) = boundaryConditions1;

%Create Matrices needed to solve
df = buildDerivative(f,stepsize);
dx = zeros(stepnumber);
for i = 1:stepnumber
    dx(i,i) = stepsize(i);
end

%Solve for f
Q = df' * dx * df;
b = zeros(1, stepnumber+1);
b(1) = boundaryConditions0;
b(length(b)) = boundaryConditions1;
b = b * Q;
Q_red = Q(2:stepnumber, 2:stepnumber);
f_unknown = Q_red \ -b(2:length(b) - 1)';
f = [boundaryConditions0; f_unknown; boundaryConditions1]

%Using eigenvalues and lagrange multipliers
[eigenvectors, eigenvalues] = eig(Q);
[m,n] = size(eigenvectors);
partition = [0];
for i = 1:size(stepsize)
    partition = [partition partition(length(partition)) + stepsize(i)];
end

%Plot eigenvectors
hold on
for i = 1:n
    plot(partition, eigenvectors(:, i), 'DisplayName', sprintf('Eigenvector of eigenvalue %f', eigenvalues(i,i) ))
end
title('Laplace Equations with ||f|| = 1')
legend
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Functions
%Input: column function f, and stepsize column vector
%Output: n-1 x n derivative matrix where n is the size of f
function derivativeMatrix = buildDerivative(f, stepsize)
    [n, ~] = size(f);
    derivativeMatrix = zeros(n-1, n);
    for i = 1:n-1
        derivativeMatrix(i,i) = -1;
        derivativeMatrix(i, i+1) = 1;
    end
    [m, ~] = size(stepsize);
    stepMatrix = zeros(m);
    for i = 1:m
        stepMatrix(i,i) = 1/stepsize(i);
    end
    derivativeMatrix = stepMatrix * derivativeMatrix;
end

%Testing again, hello world