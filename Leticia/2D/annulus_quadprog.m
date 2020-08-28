% annulus with boundary condition u = x at inner boundary and u = y at outer boundary

% create annulus
[V,F,bo_in,bo_out] = annulus(100,2);

% vectorize matrix F to obtain row location of vertices in matrix V
S = unique(reshape(F,[],1));

% write dΩ the vector of vertice indices on the boundary
bo = [bo_in; bo_out];

% find Ω the vector of vertice indices on the interior
bi = setdiff(S,bo);

% write g the vector of boundary condition(s) 
g_in = V(bo_in,1);
g_out = V(bo_out,2);
g = [g_in; g_out];

% write ui the solution vector on the interior
u = zeros(length(S),1);
u(bo) = g(S(1:size(g)));

% cotagent matrix using gptoolbox
L = cotmatrix(V,F);

% extract matrix LΩΩ from matrix L i.e. the part of L acting on bi
Lbi = -L(bi,bi);

% find the term LΩdΩ*g := b
%
% multiply with 0 entries replacing bi
b = L*u;
% then select only the corresponding entries to bi
b = b(bi);

% solve using quadprog
ubi = quadprog(Lbi,b);

% insert solution values on the interior into solution vector
u(bi) = ubi(S(1:size(ubi)));

% plot solution
plot = tsurf(F,[V u],fphong, falpha(1,0));

colormap(cbrewer('RdYlBu',40));
colorbar();

title('\nabla^2u = 0, with u = x on inner boundary, u = y on outer boundary');

axis equal
grid off;
axis off;
camlight;
add_isolines(plot);

