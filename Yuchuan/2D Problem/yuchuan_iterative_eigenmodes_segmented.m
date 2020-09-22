%% Computes iterative eigenmodes using segmented representation



%% Create mesh
[x, y] = meshgrid(1:10, 1:10);
points = [x(:), y(:)];
V = [points zeros(size(points,1),1)];
F = delaunay(points);

%% Initialization 
GMG = face_GMG(V,F);
M = face_area_matrix(V,F);
D = discontinuity(V,F);

%% Iterative method (with continuity constraint)

% n_iter =  3; 
% U = zeros(3*size(F,1),3*size(F,1)); %i_th col stores the i_th eigenmode
% 
% c = sparse(1,1,1,3*size(F,1),1); %OR
% c = c/sqrt(((c')*M*c));
% 
% options = optimset('MaxIter',200,'TolX',1e-8);
% for i = 1:n_iter
%     
%     u_old = sparse(3*size(F,1),1);
%     u_new = sparse(1,1,1,3*size(F,1),1);
%     
%     while (u_new - u_old)' * M * (u_new - u_old) > 1e-4
%         u_old = u_new;
%         Aeq = [(c')*M ; (U')*M ; D];
%         beq = sparse(1,1,1,3*size(F,1)+1 + size(D,1),1);
%         u = quadprog(GMG,[],[],[],Aeq,beq,[],[],[],options);
%         u_new = u/sqrt(((u')*M*u));
%         c = u_new;
%     end
%     
%     U(:,i) = u_new;
%     ker = null((U')*M);
%     c = ker(:,1);
%     c = c/sqrt(((c')*M*c));
%     
% end



%% Iterative method (without continuity constraint)

n_iter =  3; 
U = zeros(3*size(F,1),3*size(F,1)); %i_th col stores the i_th eigenmode

c = sparse(1,1,1,3*size(F,1),1); %OR
c = c/sqrt(((c')*M*c));

options = optimset('MaxIter',200,'TolX',1e-8);
for i = 1:n_iter
    
    u_old = sparse(3*size(F,1),1);
    u_new = sparse(1,1,1,3*size(F,1),1);
    it = 0;
    
    while (u_new - u_old)' * M * (u_new - u_old) > 1e-8 && it <= 400
        u_old = u_new;
        Aeq = [(c')*M ; (U')*M];
        beq = sparse(1,1,1,3*size(F,1)+1,1);
        u = quadprog(GMG,[],[],[],Aeq,beq,[],[],[],options);
        u_new = u/sqrt(((u')*M*u));
        c = u_new;
        it = it + 1;
    end
    
    U(:,i) = u_new;
    ker = null((U')*M);
    c = ker(:,1);
    c = c/sqrt(((c')*M*c));
    
end

%% Plot eigenmodes
face_plotting(V,F,U(:,1));



