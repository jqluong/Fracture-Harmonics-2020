function u = laplace_eq_2D(V, F, B)
    % solve laplace's equation on a 2D tiangle mesh surface
    %
    % V is vertex list
    %
    % F is face list
    %
    % B is boundary conditions. B should be a k by 2 matrix, where the
    % first column is vertex position in V, and
    % second column is the value of u at the corresponding (x, y)
    %
    % output is the solution u of \Delta u = 0 given boundary condition B
    
    % generate cotangent matrix
    % notation: L = cotangent matrix, subscript b for boundary (corresponds
    % to \partial in notes), and subscript i for interior (corresponds to
    % \Omega in notes)
    L = cotmatrix(V,F);
    % find column indices corresponding to non-boundary condition vertices
    nbv = setdiff(1:length(L),B(:,1));
    % find submatrix that operates only on non-boundary condition vertices
    L_ii = L(nbv, nbv);
    
    % set up u with boundary condition
    u = zeros(length(V),1);
    u(B(:,1)) = B(:,2);

    % if u_int = 0, L*u will have contributions from L_ib*u_b and L_bb*u_b,
    % so perform multiplication and extract values corresponding to 
    % non-boundary condition points (L_ib*u_b)
    b = L*u;
    b = b(nbv);
    
    u_int = L_ii \ b;
    
    u(nbv) = -u_int;
end
