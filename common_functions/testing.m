% script to test various utilities for segment based functions
x = [0;1;2;3;1;2;3;4];
y = [2;1;1;1;0;3;1;0];
n = length(y)/2; % 4

G = segment_derivative(n);
D = segment_discontinuity(n);
M = segment_midpoint(n);

disp('derivative');
G*y

disp('discontinuity');
D*y

disp('midpoint');
M*y

disp('l1 norm');
segment_norm(y,1,x)
disp('l2 norm');
segment_norm(y,2,x)

hold on
segment_plot(x,y)
segment_plot(x,y+1,'r')

% note curly brace for string array, square for color
segment_plot_legend({'f1';'f2'}, ['b';'r']); % note curly brace for string array, square for color
hold off