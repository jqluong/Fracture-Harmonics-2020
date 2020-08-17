%% script to demonstrate various utilities for segment based functions %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function
x = [0;1;2;3;1;2;3;4];
y = [2;1;1;1;0;3;1;0];
n = length(y)/2; % 4

%% matrix generating functions - derivative, discontinuity, midpoint
G = segment_derivative(n);
D = segment_discontinuity(n);
M = segment_midpoint(n);

disp('derivative');
G*y

disp('discontinuity');
D*y

disp('midpoint');
M*y

%% function norms
disp('l1 norm');
segment_norm(y,1,x)
disp('l2 norm');
segment_norm(y,2,x)
segment_norm(y,2)% if x is even spacing on [0,1]

disp('segment metric');
segment_metric(x,y,y.^2 - y + 1)

%% plotting
hold on
segment_plot(x,y) % color defaults to blue
segment_plot(x,y+1,'r')
segment_plot(x, (y.^2+1)-y, 'c')

% note curly brace for string array, square for color
segment_plot_legend({'f';'f + 1';'f^2 - f + 1'}, ['b';'r';'c']);
title('testing of segment\_plot and segment\_plot\_legend');
hold off