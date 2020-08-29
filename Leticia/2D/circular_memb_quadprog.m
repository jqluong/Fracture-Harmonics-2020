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
g = zeros(length(B));


% call function to solve using quadprog
u = laplace_eq_2D_quadprog(V,F,[B,g]);

% plot solution
plot_surf(F,V,u);