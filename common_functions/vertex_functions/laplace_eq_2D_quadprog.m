function u = laplace_eq_2D_quadprog(V, F, B)
    % solve laplace's equation on a 2D triangle mesh surface using quadprog
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
    
    % find number of vertices
    [n,~] = size(V);
    
    
    % find number of boundary vertices
    [k,~] = size(B);
    
    
    % generate cotagent matrx
    L = - cotmatrix(V,F);
    
    % generate boundary condition vector
    g(1:k) = B(1:k,2);
    
    % matrix used in equality constraint Aeq*u = g
    
    Aeq_i = zeros(k,1);
    Aeq_j = zeros(k,1);
    Aeq_v = zeros(k,1);
    
    for m = 1:k
        Aeq_i(m) = m;
        Aeq_j(m) = B(m,1);
        Aeq_v(m) = 1;
    end
    
    Aeq = sparse(Aeq_i,Aeq_j,Aeq_v,k,n);
    
    % solve using quadprog
    u = quadprog(L, zeros(n,1), [],[],Aeq,g);
    
end