opengl software
clear all
n = input('Number of Points: ');
x = 0:0.1:1;
y = x';
f = x + y;
[dx, dy] = gradient(f,0.1);
gradient = dx.^2 + dy.^2;
dirichletEnergy = (1/n)^2*sum(gradient, 'all')
