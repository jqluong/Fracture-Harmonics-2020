% annulus with boundary condition u = x at inner boundary and u = y at
% outer boundary using quadprog


% generate annulus
[V,F,bo_in,bo_out] = annulus(100,2);

[m,~] = size(bo_in);

[n,~] = size(bo_out);


% generate boundary conditions matrix
B = zeros(m+n,2);

% write B the matrix of boundary condition(s) 
% indice values on the first column
B(1:m,1) = bo_in;
B(m+1:m+n,1) = bo_out;
% function values on the second column
B(1:m,2) = V(bo_in,1);        % u = x at inner boundary
B(m+1:m+n,2) = V(bo_out,2);      % u = y at outer boundary

% call function to solve using quadprog
u = laplace_eq_2D_quadprog(V,F,B);

% plot solution
plot_surf(F,V,u);