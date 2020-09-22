%% Computes iterative eigenmodes using segmented representation



%% Create mesh
[x, y] = meshgrid(1:10, 1:10);
points = [x(:), y(:)];
V = [points zeros(size(points,1),1)];
F = delaunay(points);

%% Initialization 
L = -1*cotmatrix(V,F);
M = massmatrix(V,F,'barycentric');


%% Iterative method
n_iter =  6; 
U = zeros(size(V,1),size(V,1)); %i_th col stores the i_th eigenmode

c = sparse(1,1,1,size(V,1),1); 
c = c/sqrt(((c')*M*c));

options = optimset('MaxIter',200,'TolX',1e-8);
for i = 1:n_iter
    
    u_old = sparse(size(V,1),1);
    u_new = ones(size(V,1),1);
    
    while (u_new - u_old)' * M * (u_new - u_old) > 1e-4
        u_old = u_new;
        Aeq = [(c')*M ; (U')*M];
        beq = sparse(1,1,1,size(V,1)+1,1);
        u = quadprog(L,[],[],[],Aeq,beq,[],[],[],options);
        u_new = u/sqrt(((u')*M*u));
        c = u_new;
    end
    
    U(:,i) = u_new;
    ker = null((U')*M);
    c = ker(:,1);
    c = c/sqrt(((c')*M*c));
    
end

%% Plot eigenmodes

X = reshape(V(:,1),10,10)';
Y = reshape(V(:,2),10,10)';

figure
subplot(2,3,1)
surf(X,Y,reshape(U(:,1),10,10)')
subplot(2,3,2)
surf(X,Y,reshape(U(:,2),10,10)')
subplot(2,3,3)
surf(X,Y,reshape(U(:,3),10,10)')
subplot(2,3,4)
surf(X,Y,reshape(U(:,4),10,10)')
subplot(2,3,5)
surf(X,Y,reshape(U(:,5),10,10)')
subplot(2,3,6)
surf(X,Y,reshape(U(:,6),10,10)')
