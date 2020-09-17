%% Computes iterative eigenmodes using segmented representation



%% Create mesh
[x, y] = meshgrid(1:10, 1:10);
points = [x(:), y(:)];
V = [points zeros(size(points,1),1)];
F = delaunay(points);

%% Initialization 
GMG = face_GMG(V,F);
M = face_area_matrix(V,F);

%% Iterative method

n_iter =  5; %size(F,1)+5;  %number of iterations
U = zeros(3*size(F,1),3*size(F,1)); %i_th col stores the i_th eigenmode


for i = 1:n_iter
    ker = null((U')*M);
    c = ker(:,1);
    Aeq = [(c')*M ; (U')*M];
    beq = sparse(1,1,1,3*size(F,1)+1,1);
    u = quadprog(GMG,[],[],[],Aeq,beq);
    u_normalized = u/sqrt(((u')*M*u));
    U(:,i) = u_normalized;
end


%% Plot eigenmodes
face_plotting(V,F,U(:,4));


%% Check energy of eigenmode
%energy = U(:,2)' * GMG * U(:,2);

