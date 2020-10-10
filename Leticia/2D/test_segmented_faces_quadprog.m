<<<<<<< Updated upstream
% generate circle
s = 100;
theta = [0:2*pi()/s:2*pi()]';
x = cos(theta);
y = sin(theta);
=======
% construct mesh
f = 10;
x = 0:1/f:1;
y = 0:1/f:1;

[X,Y] = meshgrid(x,y);
>>>>>>> Stashed changes

V = [X(:), Y(:)];

F = delaunay(V);

[m,~] = size(V);

V_Z = zeros(m,1);


u = laplace_eq_2D_seg_quadprog([ V V_Z ],F);

face_plot(V,F,u)