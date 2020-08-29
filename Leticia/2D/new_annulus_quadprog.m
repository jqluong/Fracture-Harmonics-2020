% annulus with boundary condition u = x at inner boundary and u = y at
% outer boundary using quadprog


% generate annulus
[V,F,bo_in,bo_out] = annulus(100,2);

[m,~] = size(bo_in);

[n,~] = size(bo_out);


% generate boundary conditions matrix
%B = zeros(m+n,2);

% write B the matrix of boundary condition(s) 
% indice values on the first column
%B(1:m,1) = bo_in;
%B(m+1:m+n,1) = bo_out;
% function values on the second column
%B(1:m,2) = V(bo_in,1);        % u = x at inner boundary
%B(m+1:m+n,2) = V(bo_out,2);      % u = y at outer boundary


% sanity check using previously perfomed boundary condition
% collect points on boundary
B = unique(reshape(outline(F),[],1));

% determine inner boundary or outer boundary
B1 = B(abs(normrow(V(B,:)) - 2) <= 0.1); % outer
B2 = B(abs(normrow(V(B,:)) - 1) <= 0.1); % inner
B2p = B2(V(B2, 2) >= 0); % inner; y >= 0
B2m = B2(V(B2, 2) < 0); % inner; y < 0
B = [B1; B2p; B2m];

g1 = sin(4*(atan(V(B1,1) ./ V(B1,2))+pi()/8));
g2p = ones(length(B2p),1);
g2m = -ones(length(B2m),1);
g = [g1; g2p; g2m];

% call function to solve using quadprog
u = laplace_eq_2D_quadprog(V,F,[B,g]);

% plot solution
plot_surf(F,V,u);

% title for sanity check
title('\nabla^2u = 0, with u = cos(4\theta) on outer boundary, u = \pm 1 on inner boundary');