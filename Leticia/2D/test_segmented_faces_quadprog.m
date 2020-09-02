
% generate circle
s = 100;
theta = [0:2*pi()/s:2*pi()]';
x = cos(theta);
y = sin(theta);

[V, F] = triangle([x,y], 'Quality', 30, 'MaxArea', 0.001);


u = laplace_eq_2D_seg_quadprog(V,F);

face_plot(V,F,u)