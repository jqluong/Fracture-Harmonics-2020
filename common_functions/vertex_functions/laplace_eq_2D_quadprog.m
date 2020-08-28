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
    L = cotmatrix(V,F);
    
    % generate boundary condition vector
    g = zeros(n,1);
    g(B(1:k,1)) = B(1:k,2);

    
    % solve using quadprog
    u = quadprog(L, zeros(n,1), [],[],speye(n),g);
    
end