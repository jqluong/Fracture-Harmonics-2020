
function [V,F] = mesh_circle(s)
% generates circular mesh
% s: number of partitions on interval from 0 to 2*pi

theta = [0:2*pi()/s:2*pi()]';
x = cos(theta);
y = sin(theta);

[V, F] = triangle([x,y], 'Quality', 30, 'MaxArea', 0.001);
end