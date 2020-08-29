
% torus with boundary condition u = x + y at z = 0

% read .obj file
[V,F] = readOBJfast('torusNF.obj');

x = V(1:n,1);
y = V(1:n,2);
z = V(1:n,3);

% vectorize matrix F to obtain row location of vertices in matrix V
S = unique(reshape(F,[],1));

% find dΩ vertice indices, recall boundary is z = 0, here R = 1.25, r =
% 0.75
bo_in = S(normrow(V) <= 0.79); 
bo_out = S(1.22 <= normrow(V));

bo = [bo_in; bo_out];


% write down boundary condition
g = V(bo,1) + V(bo,2);
%g = zeros(length(bo),1);


% cotagent matrix using gptoolbox
L = cotmatrix(V,F);


% call quadprog to find solution (ignore '2D' in the function name)
u = laplace_eq_2D_quadprog(V,F,[bo,g]);

% plot solution

%scatter3(x,y,z,40,u,'filled')    % draw the scatter plot
%ax = gca;
%ax.XDir = 'reverse';
%view(-2,2)
%xlabel('x')
%ylabel('y')
%zlabel('z')
%axis equal

%cb = colorbar;  

p = patch('Faces',F,'Vertices',V,'FaceVertexCData',u);
p.FaceColor = 'interp';
colorbar
axis equal
title('\nabla^2u = 0, with u = x + y at z = 0');





