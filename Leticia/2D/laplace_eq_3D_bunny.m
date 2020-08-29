
% bunny with boundary condition u = 0 'around' (x,y,z) =
% (0.08513,0.8906,-0.3618), (x,y,z) = (-0.6699,0.8745,-0.7278) and (x,y,z)
% = (0.9831,-0.5556,0.2712)

% read .obj file
[V,F] = readOBJfast('bunny.obj');

x = V(1:n,1);
y = V(1:n,2);
z = V(1:n,3);

% vectorize matrix F to obtain row location of vertices in matrix V
S = unique(reshape(F,[],1));

% find dΩ vertice indices
bo = S(normrow(V) <= 0.97 & 0.96 <= normrow(V) & V(:,2) <= 0.9 & 0.8 <= V(:,2)); 
bo2 = S(normrow(V) <= 1.33 & 1.31 <= normrow(V) & V(:,2) <= 0.9 & 0.85 <= V(:,2)); 
bo3 = S(normrow(V) <= 1.17 & 1.15 <= normrow(V) & V(:,2) <= -0.45 & -0.65 <= V(:,2)); 



bo = [bo1;bo2;bo3];


% write down boundary condition
g = zeros(length(bo),1);

% cotagent matrix using gptoolbox
L = cotmatrix(V,F);


% call quadprog to find solution (ignore '2D' in the function name)
u = laplace_eq_2D_quadprog(V,F,[bo,g]);

% plot solution
p = patch('Faces',F,'Vertices',V,'FaceVertexCData',u);
p.FaceColor = 'interp';
colorbar
axis equal
title('\nabla^2u = 0, with u = 0 at the tip of ears and tail');
