% circular membrane with boundary condition u = 0 using quadprog

% generate circle
s = 100;
theta = [0:2*pi()/s:2*pi()]';
x = cos(theta);
y = sin(theta);

[V, F] = triangle([x,y], 'Quality', 30, 'MaxArea', 0.001);


% write B the matrix of boundary condition(s) 
B = unique(reshape(outline(F),[],1));

% create vector with function values at boundary
%g = zeros(length(B));

% sanity check using previously perfomed boundary condition
BP1 = B(V(B,2) <= 0);
BP2 = B(V(B,2) > 0);
B = [BP1; BP2];

g1 = V(BP1, 1).^3;
g2 = V(BP2, 2).^2;
g = [g1; g2];


% call function to solve using quadprog
u = laplace_eq_2D_quadprog(V,F,[B,g]);

% plot solution
plot_surf(F,V,u);

% title for sanity check
title('\nabla^2u = 0, with u = x^3 at y < 0 and u = y^2 at y > 0 on unit circle ');