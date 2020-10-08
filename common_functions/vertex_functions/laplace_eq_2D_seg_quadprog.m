function u = laplace_eq_2D_seg_quadprog(V, F)
    % solve laplace's equation on a 2D triangle mesh surface using quadprog
    % for 'segmented' functions, i.e. functions with face edges that may
    % not match
    %
    % V is matrix of vertex position
    %
    % F is matrix of face indices
    %
    % B is boundary conditions. B should be a k by 2 matrix, where the
    % first column is vertex position in V, and
    % second column is the value of u at the corresponding (x, y)
    %
    % output is the solution u of \Delta u = 0 given boundary condition B
    %
    % the format of solution u is u = [ u(F1_v1) u(F1_v2) u(F1_v3) u(F2_v1) u(F2_v2) u(F2_v3) ... ]
    
    E = edges(F); % edges matrix
    
    [m,~] = size(E); 
    
    [k,~] = size(F);
   
    % generate gradient matrix to act on [ u t ]
    cellArea = repelem(power(doublearea(V, F)/2, 0.5), 2);
    %M = spdiags(cellArea, 0, k*2, k*2);
    
    G = transpose(face_grad(V,F)) * speye(k*2,k*2) * face_grad(V,F);  % size(G) = 3k x 3k
    tf = issymmetric(G)
    O_right = zeros(3*k,m); 
    O_bottom = zeros(m,3*k+m);
    G_tilde = [ G O_right; O_bottom ];
    
    % generate discontinuity matrix
    V = [V zeros(length(V))];
    D = discontinuity(V,F);
    V = V(:,1:2);
    
    % generate vector f used in the constraint
    u_placeholder = zeros(3*k,1);
    t_placeholder = ones(m,1);
    f = [ u_placeholder; t_placeholder ];
    
    % generate |E|x|E| identity matrix used in the constraint block matrix
    I = speye(m);
    
    % generate inequality constraint matrix and vector
    A = [ -D -I; D -I];
    b = zeros(2*m,1);
    
    %boundary conditions
    B = unique(reshape(outline(F),[],1));
    BP1 = B(V(B,1) == 0);
    BP2 = B(V(B,1) == 1);
    
    %bad way to get locations of boundary vertices in face based function
    %synatx fix later
    BP1new = [];
    for i = 1:length(BP1)
        [r,s] = size(F);
        for j = 1:r
            for k = 1:s
                if F(j,k) == BP1(i)
                    BP1new = [BP1new 3*(j-1) + k];
                end
            end
        end
        
    end
    BP1new = BP1new';
    BP2new = [];
    for i = 1:length(BP2)
        [r,s] = size(F);
        for j = 1:r
            for k = 1:s
                if F(j,k) == BP2(i)
                    BP2new = [BP2new 3*(j-1) + k];
                end
            end
        end
        
    end
    BP2new = BP2new';
    [k,~] = size(F);
    
    equalityConstraint = sparse(BP2new,ones(length(BP2new),1),ones(length(BP2new),1), 3*k + m, 1);
    equalityMatrix = sparse([BP1new; BP2new], [BP1new; BP2new], ones(length(BP1new) + length(BP2new), 1), 3*k + m, 3*k + m);
    
    y = quadprog(G_tilde,f,A,b,equalityMatrix,equalityConstraint);
    
    u = y(1:3*k);
    %Constructing a function for debugging
    v = ones(length(u),1);
    for i = 1:length(BP1)
        [r,~] = size(F);
        for j = 1:r
            if (BP1(i) == F(j,1) || BP1(i) == F(j,2) || BP1(i) == F(j,3))
                v(3*(j-1) + 1) = 0;
                v(3*(j-1) + 2) = 0;
                v(3*(j-1) + 3) = 0;
            end
        end
    end
norm(face_grad(V,F)*v,2) + norm(D*v,1)
norm(face_grad(V,F)*u,2) + norm(D*u,1)
end