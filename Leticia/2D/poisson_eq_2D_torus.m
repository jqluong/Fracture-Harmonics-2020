
% torus with boundary condition u = x + y at z = 0

% read .obj file
[V,F] = readOBJfast('torusNF.obj');

% vectorize matrix F to obtain row location of vertices in matrix V
S = unique(reshape(F,[],1));

% find dΩ vertice indices, recall boundary is z = 0
bo = S(V(S,3) == 0); 


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





