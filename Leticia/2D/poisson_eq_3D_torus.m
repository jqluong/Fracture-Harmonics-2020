
% torus with boundary condition u = x + y at z = 0

% read .obj file
[V,F] = readOBJfast('torusNF.obj');

tsurf(F,V)

% vectorize matrix F to obtain row location of vertices in matrix V
S = unique(reshape(F,[],1));

% find dΩ vertice indices, recall boundary is z = 0, here R = 1.25, r =
% 0.75
bo_in = S(normrow(V) <= 0.79); 
bo_out = S(1.22 <= normrow(V));

bo = [bo_in; bo_out];

% find Ω vertice indices
bi = setdiff(S,bo);


% write down boundary condition
g = V(bo,1) + V(bo,2);

% write down vector of solution
u = zeros(length(S),1);
u(bo) = g(S(1:size(g)));

% cotagent matrix using gptoolbox
L = cotmatrix(V,F);

% extract from matrix LΩΩ i.e. the part of L acting on bi
Lbi = L(bi,bi);

% find the term LΩdΩ*g := b
%
% multiply with 0 entries replacing bi
b = L*u;
% then select only the corresponding entries to bi
b = b(bi);

% solve
ubi = Lbi \ b;

u(bi) = -ubi;





