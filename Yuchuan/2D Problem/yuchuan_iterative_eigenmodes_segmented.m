%% Computes iterative eigenmodes using segmented representation



%% Create mesh
[x, y] = meshgrid(1:10, 1:10);
points = [x(:), y(:)];
V = [points zeros(size(points,1),1)];
F = delaunay(points);



%% Initialization 
GMG = face_GMG(V,F);
M = face_massmatrix(V,F);
D = discontinuity(V,F);
area_matrix = face_area_matrix(V,F);










%% (1) argmin u^t * GMG * u (with continuity constraint)
% n_iter =  2; 
% U = zeros(3*size(F,1),3*size(F,1)); %i_th col stores the i_th eigenmode
% c = sparse(1,1,1,3*size(F,1),1); %OR
% c = c/sqrt(((c')*M*c));
% options = optimset('MaxIter',200,'TolX',1e-8);
% for i = 1:n_iter
%     u_old = sparse(3*size(F,1),1);
%     u_new = sparse(1,1,1,3*size(F,1),1);
%     while norm(u_new - u_old) > 1e-8
%         u_old = u_new;
%         Aeq = [(c')*M ; (U')*M ; D];
%         beq = sparse(1,1,1,3*size(F,1)+1 + size(D,1),1);
%         u = quadprog(GMG,[],[],[],Aeq,beq,[],[],[],options);
%         u_new = u/sqrt(((u')*M*u));
%         c = u_new;
%     end
%     U(:,i) = u_new;
%     ker = null((U')*M);
%     c = ker(:,1);
%     c = c/sqrt(((c')*M*c));
% end








%% (2) argmin u^t * GMG * u (without continuity constraint)
% n_iter =  3; 
% U = zeros(3*size(F,1),3*size(F,1)); %i_th col stores the i_th eigenmode
% c = sparse(1,1,1,3*size(F,1),1); %OR
% c = c/sqrt(((c')*M*c));
% options = optimset('MaxIter',200,'TolX',1e-8);
% for i = 1:n_iter
%     u_old = sparse(3*size(F,1),1);
%     u_new = sparse(1,1,1,3*size(F,1),1);
%     it = 0;
%     while norm(u_new - u_old) > 1e-8 && it <= 400
%         u_old = u_new;
%         Aeq = [(c')*M ; (U')*M];
%         beq = sparse(1,1,1,3*size(F,1)+1,1);
%         u = quadprog(GMG,[],[],[],Aeq,beq,[],[],[],options);
%         u_new = u/sqrt(((u')*M*u));
%         c = u_new;
%         it = it + 1;
%     end
%     U(:,i) = u_new;
%     ker = null((U')*M);
%     c = ker(:,1);
%     c = c/sqrt(((c')*M*c));
% end






%% (3) argmin u^t * GMG * u + mu*|D*u| 
n_iter =  4; 
U = zeros(3*size(F,1),3*size(F,1)); %i_th col stores the i_th eigenmode
U(1:6 , 1) = 0.2; %remove if not needed   %forces first eigenmode to be discontinuous
%c = sparse(1,1,1,3*size(F,1),1); %OR
c = rand(3*size(F,1),1);
c = c/sqrt(((c')*area_matrix*c));
k = size(D,1);
mu = 1e-7;
H = [GMG sparse(size(GMG,1),k) ; sparse(k,size(GMG,1)+k)];
f = [zeros(1,3*size(F,1)) ones(1,k)]' ;
f = sparse(f);
A = mu * [D (-1*speye(k)) ; -1*D (-1*speye(k))];
b = sparse(2*k,1);
beq = sparse(1,1,1,3*size(F,1)+1,1);
options = optimset('MaxIter',400,'TolX',1e-8);

for i = 2:n_iter  % change back to 1:n_iter!!!
    u_old = sparse(3*size(F,1),1);
    u_new = sparse(1,1,1,3*size(F,1),1);
    it = 0;
    while norm(u_new - u_old) > 1e-8 && it <= 400
        u_old = u_new;
        Aeq = [(c')*area_matrix sparse(1,k) ; (U')*area_matrix sparse(3*size(F,1),k)];
        u_t = quadprog(H,f,A,b,Aeq,beq,[],[],[],options);
        u = u_t(1:3*size(F,1));
        u_new = u/sqrt(((u')*area_matrix*u));
        c = u_new;
        it = it + 1;
    end
    U(:,i) = u_new;
    ker = null((U')*area_matrix);
    c = ker(:,1);
    c = c/sqrt(((c')*area_matrix*c));
end






%% Plot eigenmodes
subplot(1,3,1)
face_plotting(V,F,U(:,1));
zlim([-0.5 0.5]);

subplot(1,3,2)
face_plotting(V,F,U(:,2));
zlim([-0.5 0.5]);

subplot(1,3,3)
face_plotting(V,F,U(:,3));
zlim([-0.5 0.5]);


%% Check energies

E = U(:,3)' *GMG*U(:,3)
Du = norm(D*U(:,3) , 1)

